#include "dmrg_vbatch.h"

#define USE_MALLOC

#ifdef USE_INTEL_MKL
#include "mkl.h"
#endif

#ifdef USE_MAGMA
#include "cuda.h"
#include "cuda_runtime.h"
#include "magma_v2.h"

   static int isfirst = 1;
   static magma_queue_t queue = 0;
   static int device = 0;
#endif

#ifdef __cplusplus
extern "C" {
#endif
extern 
void dgemm_( const char *transa, const char *transb,
             const int *m, const int *n, const int *k,
             const double *alpha,
             const double *A, const int *lda,
             const double *B, const int *ldb,
             const double *beta,
                   double *C, const int *ldc );

extern
void dcopy_( const int *n, 
             const double *x, const int *incx,
                   double *y, const int *incy );

extern
void dlacpy_( const char *uplo,
              const int *m, const int *n,
              const double *A, const int *lda,
                    double *B, const int *ldb );
#ifdef __cplusplus
}
#endif

void dmrg_init()
{
#ifdef USE_MAGMA
  if (isfirst) {
     isfirst = 0;

     device = 0;
     magma_setdevice( device );
     magma_queue_create( device, &queue );
     assert( queue != 0 );
     };
#endif

}

void dmrg_dgetvector( const int n, 
                      const double *dx_src, const int incx, 
                            double *hy_dst, const int incy )
{
#ifdef USE_MAGMA
 magma_dgetvector(n, 
                  dx_src, incx, 
                  hy_dst, incy, 
                  queue );
#else
 dcopy_( &n, 
         dx_src, &incx, 
         hy_dst, &incy );
#endif
}


void dmrg_dsetvector(  const int n,  
                       const double *hx_src,  const int incx, 
                             double *dy_dst,  const int incy )
{
#ifdef USE_MAGMA
 magma_dsetvector(  n,
                    hx_src, incx,
                    dy_dst, incy,
                    queue );
#else
 dcopy_( &n, 
          hx_src, &incx, 
          dy_dst, &incy );
#endif
}


void dmrg_dgetmatrix( const int m, const int n, 
                      const double *dA_src, const int ldda, 
                      double *hB_dst, const int ldb )
{
#ifdef USE_MAGMA
  magma_dgetmatrix( m,n,
                    dA_src, ldda,
                    hB_dst, ldb,
                    queue );
#else
 const char *uplo = "A";
 dlacpy_( uplo, &m, &n, 
          dA_src, &ldda, 
          hB_dst, &ldb );
#endif
}


void dmrg_dsetmatrix( const int m, const int n, 
                      const double *hA_src, const int lda, 
                      double *dB_dst, const int lddb )
{
#ifdef USE_MAGMA
 magma_dsetmatrix( m,n,
                   hA_src, lda,
                   dB_dst, lddb,
                   queue );
#else
 const char *uplo = "A";
 dlacpy_( uplo, &m, &n, hA_src, &lda, dB_dst, &lddb );
#endif
}

void dmrg_dgemm_vbatch( char ctransa_array[],
                        char ctransb_array[],
                        int  m_array[],
                        int  n_array[],
                        int  k_array[],
                        double alpha_array[],
                        double *a_array[], int    lda_array[],
                        double *b_array[], int    ldb_array[],
                        double beta_array[],
                        double *c_array[], int    ldc_array[],
                        int group_count, 
                        int group_size[] )
{

 const int idebug = 1;


#ifdef USE_INTEL_MKL
  {

  MKL_INT m_array_mkl[group_count];
  MKL_INT n_array_mkl[group_count];
  MKL_INT k_array_mkl[group_count];

  MKL_INT lda_array_mkl[group_count];
  MKL_INT ldb_array_mkl[group_count];
  MKL_INT ldc_array_mkl[group_count];

  MKL_INT group_size_mkl[group_count];

  CBLAS_TRANSPOSE transa_array[group_count];
  CBLAS_TRANSPOSE transb_array[group_count];

  MKL_INT group_count_mkl = group_count;

  int igroup = 0;
  for(igroup=0; igroup < group_count; igroup++) {
      char transa = ctransa_array[igroup];
      char transb = ctransb_array[igroup];
      int is_transa = (transa == 'T') || (transa == 't');
      int is_transb = (transb == 'T') || (transb == 't');
  
      transa_array[igroup] = is_transa ?  CblasTrans : CblasNoTrans;
      transb_array[igroup] = is_transb ?  CblasTrans : CblasNoTrans;
 
      m_array_mkl[igroup] = m_array[igroup];
      n_array_mkl[igroup] = n_array[igroup];
      k_array_mkl[igroup] = k_array[igroup];

      lda_array_mkl[igroup] = lda_array[igroup];
      ldb_array_mkl[igroup] = ldb_array[igroup];
      ldc_array_mkl[igroup] = ldc_array[igroup];

      group_size_mkl[igroup] = group_size[igroup];
      };

  cblas_dgemm_batch( CblasColMajor, 
                     transa_array, transb_array,
                     m_array_mkl, n_array_mkl, k_array_mkl,
                     alpha_array,  (const double **) &(a_array[0]), lda_array_mkl, 
                                   (const double **) &(b_array[0]), ldb_array_mkl,
                     beta_array,   c_array, ldc_array_mkl,
                     group_count_mkl, group_size_mkl );
   };
#else
   /*
    ---------------------
    expand out the groups
    ---------------------
    */
  int batch_size = 0;
  int igroup = 0;
  for(igroup=0; igroup < group_count; igroup++) {
     batch_size += group_size[igroup];
     };

  /*
   ----------------------------------------------------------
   Note magma_dgemm_vbatched need arrays of size batch_size+1
   Set vbatch_dim to be multiple of 32 for better alignment in memory
   ----------------------------------------------------------
   */
  const int ialign = 32;
  int vbatch_dim = ialign * ICEIL( (batch_size+1),ialign );
#ifdef USE_MALLOC
  double *alpha_vbatch = (double *) dmrg_malloc( sizeof(double) * (vbatch_dim));

  double *beta_vbatch = (double *) dmrg_malloc( sizeof(double) * (vbatch_dim));
  int *m_vbatch = (int *) dmrg_malloc(sizeof(int) * (vbatch_dim));
  int *n_vbatch = (int *) dmrg_malloc(sizeof(int) * (vbatch_dim));
  int *k_vbatch = (int *) dmrg_malloc(sizeof(int) * (vbatch_dim));
  char *transa_vbatch = (char *) dmrg_malloc(sizeof(char) *(vbatch_dim));
  char *transb_vbatch = (char *) dmrg_malloc(sizeof(char) *(vbatch_dim));
  int *lda_vbatch = (int *) dmrg_malloc(sizeof(int)*(vbatch_dim));
  int *ldb_vbatch = (int *) dmrg_malloc(sizeof(int)*(vbatch_dim));
  int *ldc_vbatch = (int *) dmrg_malloc(sizeof(int)*(vbatch_dim));

  double **a_vbatch = (double **) dmrg_malloc( sizeof(double *) * (vbatch_dim));
  double **b_vbatch = (double **) dmrg_malloc( sizeof(double *) * (vbatch_dim));
  double **c_vbatch = (double **) dmrg_malloc( sizeof(double *) * (vbatch_dim));

  assert( alpha_vbatch != NULL );
  assert( beta_vbatch != NULL );

  assert( m_vbatch != NULL );
  assert( n_vbatch != NULL);
  assert( k_vbatch != NULL);

  assert( transa_vbatch != NULL );
  assert( transb_vbatch != NULL );

  assert( lda_vbatch != NULL );
  assert( ldb_vbatch != NULL );
  assert( ldc_vbatch != NULL );

  assert( a_vbatch != NULL );
  assert( b_vbatch != NULL );
  assert( c_vbatch != NULL );

#else
  double alpha_vbatch[vbatch_dim];
  double beta_vbatch[vbatch_dim];
  int m_vbatch[vbatch_dim];
  int n_vbatch[vbatch_dim];
  int k_vbatch[vbatch_dim];
  char transa_vbatch[vbatch_dim];
  char transb_vbatch[vbatch_dim];
  int lda_vbatch[vbatch_dim];
  int ldb_vbatch[vbatch_dim];
  int ldc_vbatch[vbatch_dim];

  double **a_vbatch = a_array;
  double **b_vbatch = b_array;
  double **c_vbatch = c_array;

#endif

  int idx = 0;
  for(igroup=0; igroup < group_count; igroup++) {
      int i = 0;
      for(i=0; i < group_size[igroup]; i++) {
         transa_vbatch[idx] =  ctransa_array[igroup];
         transb_vbatch[idx] =  ctransb_array[igroup];

         m_vbatch[idx] = m_array[igroup];
         n_vbatch[idx] = n_array[igroup];
         k_vbatch[idx] = k_array[igroup];

         alpha_vbatch[idx] = alpha_array[igroup];
         beta_vbatch[idx] = beta_array[igroup];
         
         lda_vbatch[idx] = lda_array[igroup];
         ldb_vbatch[idx] = ldb_array[igroup];
         ldc_vbatch[idx] = ldc_array[igroup];

         a_vbatch[idx] = a_array[idx];
         b_vbatch[idx] = b_array[idx];
         c_vbatch[idx] = c_array[idx];

         idx = idx + 1;
         };
      };
      

#ifdef USE_MAGMA
   {
   /*
    ---------------------------------------------------------------
    Note magma_dgemm_vbatched assumes all transa, transb, alpha, beta are the same
    ---------------------------------------------------------------
    */
   int is_same_transa = 0;
   int is_same_transb = 0;
   int is_same_alpha = 0;
   int is_same_beta = 0;

   for(idx=0; idx < batch_size; idx++) {
      is_same_transa = (transa_vbatch[0] == transa_vbatch[idx]);
      if (!is_same_transa) { break; };
      };
   for(idx=0; idx < batch_size; idx++) {
      is_same_transb = (transb_vbatch[0] == transb_vbatch[idx]);
      if (!is_same_transb) { break; };
      };
   for(idx=0; idx < batch_size; idx++) {
      is_same_alpha = (alpha_vbatch[0] == alpha_vbatch[idx]);
      if (!is_same_alpha) { break; };
      };
   for(idx=0; idx < batch_size; idx++) {
      is_same_beta = (beta_vbatch[0] == beta_vbatch[idx]);
      if (!is_same_beta) { break; };
      };

   int is_same_all = is_same_transa &&
                     is_same_transb &&
                     is_same_alpha &&
                     is_same_beta;

   if (!is_same_all) {
    if (!is_same_transa) {
       fprintf(stderr,"dmrg_vbatch: is_same_transa = %d\n",
                                    is_same_transa );
       };
    if (!is_same_transb) {
       fprintf(stderr,"dmrg_vbatch: is_same_transb = %d\n",
                                    is_same_transb );
       };
    if (!is_same_alpha) {
       fprintf(stderr,"dmrg_vbatch: is_same_alpha = %d\n",
                                    is_same_alpha );
       };
    if (!is_same_beta) {
       fprintf(stderr,"dmrg_vbatch: is_same_beta = %d\n",
                                    is_same_beta );
       };
    };
   assert( is_same_all );


   int isTransA = ((transa_vbatch[0] == 'T') || (transa_vbatch[0] == 't'));
   int isTransB = ((transb_vbatch[0] == 'T') || (transb_vbatch[0] == 't'));

   magma_trans_t transA = isTransA ? MagmaTrans : MagmaNoTrans;
   magma_trans_t transB = isTransB ? MagmaTrans : MagmaNoTrans;
   
   double alpha = alpha_vbatch[0];
   double beta = beta_vbatch[0];


     

   magmablas_dgemm_vbatched( transA, transB,
           m_vbatch, n_vbatch, k_vbatch,
           alpha,
           a_vbatch, lda_vbatch, b_vbatch, ldb_vbatch,
           beta,
           c_vbatch, ldc_vbatch,
           batch_size, queue );



  };
#else
  {
  int i = 0;
  #pragma omp parallel for private(i) schedule(dynamic)
  for(i=0; i < batch_size; i++) {
     const int mm = m_vbatch[i];
     const int nn = n_vbatch[i];
     const int kk = k_vbatch[i];
     const char transa = transa_vbatch[i];
     const char transb = transb_vbatch[i];
     const double alpha = alpha_vbatch[i];
     const double beta = beta_vbatch[i];
     const int lda = lda_vbatch[i];
     const int ldb = ldb_vbatch[i];
     const int ldc = ldc_vbatch[i];
     const double *pa = a_vbatch[i];
     const double *pb = b_vbatch[i];
           double *pc = c_vbatch[i];
 
     dgemm_( &transa, &transb, 
             &mm, &nn, &kk, 
             &alpha, pa, &lda, pb, &ldb,
             &beta,  pc, &ldc );
     };
   };
#endif

#ifdef USE_MALLOC
    

  dmrg_free( alpha_vbatch);
  dmrg_free( beta_vbatch);

  dmrg_free(m_vbatch);
  dmrg_free(n_vbatch);
  dmrg_free(k_vbatch);

  dmrg_free( transa_vbatch );
  dmrg_free( transb_vbatch );

  dmrg_free( lda_vbatch );
  dmrg_free( ldb_vbatch );
  dmrg_free( ldc_vbatch );

  dmrg_free( a_vbatch );
  dmrg_free( b_vbatch );
  dmrg_free( c_vbatch );

#endif

#endif
}

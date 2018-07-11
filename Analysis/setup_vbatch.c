#include "dmrg_vbatch.h"
#include "test_vbatch.h"
void setup_vbatch( int noperator,
                   int npatches,
                   int left_patch_size_[],
                   int right_patch_size_[],
                   long *pleft_patch_start[],
                   long *pright_patch_start[],
                   long *pxy_patch_start[],
                   double *pAbatch[], int *pld_Abatch,
                   double *pBbatch[], int *pld_Bbatch,
                   int interaction_matrix_[] )
#define left_patch_size(ipatch) left_patch_size_[ (ipatch)-1]
#define right_patch_size(ipatch) right_patch_size_[ (ipatch)-1]
#define left_patch_start(ipatch) left_patch_start_[ (ipatch)-1]
#define right_patch_start(ipatch) right_patch_start_[ (ipatch)-1]
#define xy_patch_start(ipatch) xy_patch_start_[ (ipatch)-1]

#define interaction_matrix(ipatch,jpatch) interaction_matrix_[indx2f(ipatch,jpatch,npatches)]

#define hAbatch(i,j)  hAbatch_[indx2f(i,j,ld_Abatch)]
#define hBbatch(i,j)  hBbatch_[indx2f(i,j,ld_Bbatch)]
{


  long left_max_state = 0;
  long right_max_state = 0;
  {
   int ipatch = 0;
   for(ipatch=1; ipatch <=  npatches; ipatch++) {
      left_max_state += left_patch_size(ipatch);
      right_max_state += right_patch_size(ipatch);
      };
  }

  int nrow_Abatch = left_max_state;
  int ncol_Abatch = left_max_state * noperator;

  int nrow_Bbatch = right_max_state;
  int ncol_Bbatch = right_max_state * noperator;

  const int ialign = 32;
  int ld_Abatch = ialign * ICEIL( nrow_Abatch, ialign );
  int ld_Bbatch = ialign * ICEIL( nrow_Bbatch, ialign );
  
  double *Abatch_ = (double *) dmrg_malloc( sizeof(double)*ld_Abatch * 
                                                           ncol_Abatch );
  double *Bbatch_ = (double *) dmrg_malloc( sizeof(double)*ld_Bbatch * 
                                                           ncol_Bbatch );
  assert( Abatch_ != 0 );
  assert( Bbatch_ != 0 );

  int ipatch = 0;
  long *left_patch_start_ = (long *) malloc(sizeof(long)*(npatches+1));
  long *right_patch_start_ = (long *) malloc(sizeof(long)*(npatches+1));
  long *xy_patch_start_ = (long *) malloc(sizeof(long)*(npatches+1));

  assert( left_patch_start_ != NULL );
  assert( right_patch_start_ != NULL );
  assert( xy_patch_start_ != NULL );

  left_patch_start(1) = 1;
  for(ipatch=1; ipatch <= npatches; ipatch++) {
    left_patch_start(ipatch+1) = left_patch_start(ipatch) + left_patch_size(ipatch);
    };

  right_patch_start(1) = 1;
  for(ipatch=1; ipatch <= npatches; ipatch++) {
    right_patch_start(ipatch+1) = right_patch_start(ipatch) + right_patch_size(ipatch);
    };
    
  xy_patch_start(1) = 1;
  for(ipatch=1; ipatch <= npatches; ipatch++) {
     int nrowX = right_patch_size(ipatch);
     int ncolX = left_patch_size(ipatch);
     xy_patch_start(ipatch+1) = xy_patch_start(ipatch) + (nrowX*ncolX);
     };

    


  /*
   ---------------------------
   set seed to be reproducible
   ---------------------------
   */
  {
  unsigned int iseed = 13;
  srand(  iseed );
  }

 /*
  -------------------------
  fill in Abatch and Bbatch
  -------------------------
  */
  {
  long i = 0;
  long j = 0;
  const long size_Abatch = ld_Abatch * left_max_state * noperator;
  const long size_Bbatch = ld_Bbatch * right_max_state * noperator;


  {
#ifdef USE_GETSET
  double hAbatch_[size_Abatch];
#else
  double *hAbatch_ = Abatch_;
#endif

  int ipatch = 0;
  int jpatch = 0;
  
  for(jpatch=1; jpatch <= npatches; jpatch++) {
  for(ipatch=1; ipatch <= npatches; ipatch++) {
    int has_work = interaction_matrix( ipatch,jpatch);

    int ia_start = left_patch_start(ipatch);
    int ia_end   = left_patch_start(ipatch+1)-1;

    int ja_start = left_patch_start(jpatch);
    int ja_end   = left_patch_start(jpatch+1)-1;


    int i = 0;
    int j = 0;
    for(j=ja_start; j <= ja_end; j++) {
    for(i=ia_start; i <= ia_end; i++) {
       hAbatch(i,j) = 0;
       if (has_work) {
          double dval = ((double) i + (j-1)*nrow_Abatch)/( ((double) nrow_Abatch) * ncol_Abatch );
          hAbatch(i,j) =  dval;
          };
       };
       };



  };
  };

#ifdef USE_GETSET
    {
    const int incx = 1;
    const int incy = 1;
    dmrg_dsetvector( size_Abatch, 
                     &(hAbatch_[0]), incx, 
                     Abatch_,        incy );
    }
#endif
  };





  {
#ifdef USE_GETSET
  double hBbatch_[ size_Bbatch ];
#else
  double *hBbatch_ = Bbatch_;
#endif
  int ipatch = 0;
  int jpatch = 0;
  for(jpatch=1; jpatch <= npatches; jpatch++) {
  for(ipatch=1; ipatch <= npatches; ipatch++) {
    int has_work = interaction_matrix( ipatch,jpatch);

    int ib_start = right_patch_start(ipatch);
    int ib_end   = right_patch_start(ipatch+1)-1;

    int jb_start = right_patch_start(jpatch);
    int jb_end = right_patch_start(jpatch+1)-1;

    int i = 0;
    int j = 0;
    
    for(j=jb_start; j <= jb_end; j++) {
    for(i=ib_start; i <= ib_end; i++) {
       hBbatch(i,j) = 0;
       if (has_work) {
          double dval = ((double) i+j)/( ((double) nrow_Bbatch)*ncol_Bbatch );
          hBbatch(i,j) = -dval;
       };
     };
     };

    };
    };
#ifdef USE_GETSET
    {
    const int incx = 1;
    const int incy = 1;
    dmrg_dsetvector(size_Bbatch, 
                    &(hBbatch_[0]), incx, 
                    Bbatch_,        incy );
    };
#endif

  };
  
 }


  



  *pAbatch = Abatch_;
  *pBbatch = Bbatch_;
  *pleft_patch_start = left_patch_start_;
  *pright_patch_start = right_patch_start_;
  *pxy_patch_start = xy_patch_start_;

  *pld_Abatch = ld_Abatch;
  *pld_Bbatch = ld_Bbatch;
}
          
                   

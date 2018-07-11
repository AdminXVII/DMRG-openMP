#include "dmrg_vbatch.h"
#ifdef _OPENMP
#include <omp.h>
#endif

void apply_Htarget_vbatch( 
                    int noperator,
                    int npatches, 
                    long left_patch_start_[],
                    long right_patch_start_[],
                    long xy_patch_start_[],
                    double Abatch_[], int ld_Abatch,
                    double Bbatch_[], int ld_Bbatch,
                    double X_[],
                    double Y_[])

#define left_patch_start(i) left_patch_start_[(i)-1]
#define right_patch_start(i) right_patch_start_[(i)-1]
#define xy_patch_start(i) xy_patch_start_[(i)-1]
#define Abatch(i,j) Abatch_[ indx2f(i,j,ld_Abatch) ]
#define Bbatch(i,j) Bbatch_[ indx2f(i,j,ld_Bbatch) ]
#define X(i) X_[(i)-1]
#define Y(i) Y_[(i)-1]
{
 const int ialign = 32;

 double gflops1 = (double) 0.0;
 double gflops2 = (double) 0.0;
 double time_1st_vbatch = (double) 0.0;
 double time_2nd_vbatch = (double) 0.0;

/*
 ------------------
 compute  Y = H * X
 ------------------
*/
 int ipatch = 0;
 int jpatch = 0;
 int left_max_states  = left_patch_start(npatches+1)-1;
 int right_max_states = right_patch_start(npatches+1)-1;


 int left_patch_size_[npatches];
 int right_patch_size_[npatches];
#define left_patch_size(i) left_patch_size_[(i)-1]
#define right_patch_size(i) right_patch_size_[(i)-1]

 for(ipatch=1; ipatch <= npatches; ipatch++) {
     int L1 = left_patch_start(ipatch);
     int L2 = left_patch_start(ipatch+1)-1;
    
     left_patch_size(ipatch) =  L2 - L1 + 1;
     };
 for(ipatch=1; ipatch <= npatches; ipatch++) {
    int  R1 = right_patch_start(ipatch);
    int  R2 = right_patch_start(ipatch+1)-1;
    right_patch_size(ipatch) = R2 - R1 + 1;
    };

 int ngroups = npatches;
 int ngroups_dim = ialign * ICEIL( ngroups, ialign );
 int batch_size = ngroups * noperator;
 int batch_size_dim = ialign * ICEIL( batch_size, ialign );

 double alpha_array_[ngroups_dim];
 double beta_array_[ngroups_dim];
 double *a_array_[batch_size_dim];
 double *b_array_[batch_size_dim];
 double *c_array_[batch_size_dim];

 int m_array_[ngroups_dim]; 
 int n_array_[ngroups_dim]; 
 int k_array_[ngroups_dim]; 
 int group_size_[ngroups_dim]; 
 int lda_array_[batch_size_dim];
 int ldb_array_[batch_size_dim];
 int ldc_array_[batch_size_dim];

 char transa_array_[ngroups_dim];
 char transb_array_[ngroups_dim];

#define transa_array(i) transa_array_[ (i)-1]
#define transb_array(i) transb_array_[ (i)-1]
#define m_array(i) m_array_[ (i)-1]
#define n_array(i) n_array_[ (i)-1]
#define k_array(i) k_array_[ (i)-1]
#define alpha_array(i) alpha_array_[ (i)-1]
#define beta_array(i) beta_array_[ (i)-1]
#define lda_array(i) lda_array_[ (i)-1]
#define ldb_array(i) ldb_array_[ (i)-1]
#define ldc_array(i) ldc_array_[ (i)-1]
#define a_array(i) a_array_[(i)-1]
#define b_array(i) b_array_[(i)-1]
#define c_array(i) c_array_[(i)-1]
#define group_size(i) group_size_[(i)-1]

 int nrowA = left_max_states;
 int ncolA = nrowA;
 int nrowB = right_max_states;
 int ncolB = nrowB;

 

 int nrowBX = nrowB;
 int ncolBX = (ncolA * noperator );
 int ld_BX = ialign * ICEIL(nrowBX,ialign);

 double *BX_ = (double *) dmrg_malloc( (sizeof(double) * ld_BX) * (ncolA * noperator) );
 assert( BX_ != NULL );

#define BX(i,j) BX_[ indx2f(i,j,ld_BX) ]


 ipatch = 1;
 int idx = 1;
 for(jpatch=1; jpatch <= npatches; jpatch++) {
    int igroup = jpatch;
    long j1 = xy_patch_start(jpatch);
    long j2 = xy_patch_start(jpatch+1)-1;
    int nrowX = right_patch_size(jpatch);
    int ncolX = left_patch_size(jpatch);
    assert( (j2-j1+1) == (nrowX * ncolX) );


    /*
     --------------------------------------
     XJ = reshape( X(j1:j2), nrowX, ncolX )
     --------------------------------------
     */
    double *XJ = &( X(j1) );
    int ld_XJ = nrowX;
    
    int R1 = right_patch_start(jpatch);
    int R2 = right_patch_start(jpatch+1)-1;
    int L1 = left_patch_start(jpatch);
    int L2 = left_patch_start(jpatch+1)-1;
    int kmax = noperator;
    int k = 0;

    /*
     -------------------------------
     independent DGEMM in same group
     -------------------------------
     */
    group_size(igroup) = kmax;
    for(k=1; k <= kmax; k++) {
        int offsetB = (k-1)*ncolB;
        int offsetBX = (k-1)*ncolA;
    
        
        /*
        ------------------------------------------------------------------------
        BX(1:nrowBX, offsetBX + (L1:L2)) = Bbatch(1:nrowBX, offsetB + (R1:R2) ) * 
                                             XJ( 1:(R2-R1+1), 1:(L2-L1+1));
        ------------------------------------------------------------------------
        */
        transa_array(igroup) = 'N';
        transb_array(igroup) = 'N';
        int mm = nrowBX;
        int nn = L2-L1+1;
        int kk = R2-R1+1;
        m_array(igroup) = mm;
        n_array(igroup) = nn;
        k_array(igroup) = kk;

        gflops1 += ((2.0*mm)*nn)*kk;

        alpha_array(igroup) = (double) 1;
        beta_array(igroup) = (double) 0;

        c_array(idx) = &(BX(1,offsetBX+L1));  
        ldc_array(igroup) = ld_BX;

        a_array(idx) = &(Bbatch(1,offsetB+R1)); 
        lda_array(igroup) = ld_Bbatch;

        b_array(idx) = XJ;  
        ldb_array(igroup) = ld_XJ;
        idx = idx + 1;
        };
   }; 
   /*
    ------------------
    first vbatch DGEMM
    ------------------
    */
   
#ifdef _OPENMP
   time_1st_vbatch = -omp_get_wtime();
#endif
   dmrg_dgemm_vbatch( transa_array_, transb_array_,
                      m_array_, n_array_, k_array_,
                      alpha_array_,  a_array_, lda_array_, b_array_, ldb_array_,
                      beta_array_,   c_array_, ldc_array_,
                      ngroups, group_size_ );
#ifdef _OPENMP
    time_1st_vbatch += omp_get_wtime();
#endif
    gflops1 = gflops1/(1000.0*1000.0*1000.0);
   


/*
 -------------------------------------------------
 perform computations with  Y += (BX)*transpose(A)
 -------------------------------------------------
*/
   for(ipatch=1; ipatch <= npatches; ipatch++) {
     int igroup = ipatch;
   
     long i1 = xy_patch_start(ipatch);
     long i2 = xy_patch_start(ipatch+1)-1;
   
     int R1 = right_patch_start(ipatch);
     int R2 = right_patch_start(ipatch+1)-1;
   
     int L1 = left_patch_start(ipatch);
     int L2 = left_patch_start(ipatch+1)-1;
   
     int isok = ((R2-R1+1) == right_patch_size(ipatch)) && 
                ((L2-L1+1) == left_patch_size(ipatch));
     assert(isok);

     double *YI = &(Y(i1));
     int nrowYI = R2-R1+1;
     int ld_YI = nrowYI;
     int ncolYI = L2-L1+1;
     assert( (i2-i1+1) == (nrowYI * ncolYI) );
   
     /*
        --------------------------------------------------------------------
        YI(1:(R2-R1+1),1:(L2-L1+1)) = BX( R1:R2,1:ncolBX) * 
                                         transpose( Abatch( L1:L2,1:ncolBX) );
        --------------------------------------------------------------------
      */
     group_size(igroup) = 1;
     transa_array(igroup) = 'N';
     transb_array(igroup) = 'T';
     int mm = nrowYI;
     int nn = ncolYI;
     int kk = ncolBX;
     m_array(igroup) = mm;
     n_array(igroup) = nn;
     k_array(igroup) = kk;
     gflops2 += ((2.0*mm)*nn)*kk;
     alpha_array(igroup) = (double) 1;
     beta_array(igroup) = (double) 0;
     a_array(igroup) =  &(BX(R1,1));
     lda_array(igroup) = ld_BX;
     b_array(igroup) = &(Abatch(L1,1));
     ldb_array(igroup) = ld_Abatch;
     c_array(igroup) = YI;
     ldc_array(igroup) = ld_YI;
     };
     ngroups = npatches;


     
   /*
    ------------------
    second vbatch DGEMM
    ------------------
    */
#ifdef _OPENMP
   time_2nd_vbatch = -omp_get_wtime();   
#endif
   dmrg_dgemm_vbatch( transa_array_, transb_array_,
                      m_array_, n_array_, k_array_,
                      alpha_array_,  a_array_, lda_array_, b_array_, ldb_array_,
                      beta_array_,   c_array_, ldc_array_,
                      ngroups, group_size_ );
#ifdef _OPENMP
   time_2nd_vbatch += omp_get_wtime();
   gflops2 = gflops2/(1000.0*1000.0*1000.0);

   printf("1st vbatch %f gflops (gflops1=%lf,time=%lf)\n", 
          gflops1/time_1st_vbatch,  gflops1, time_1st_vbatch );
   printf("2nd vbatch %f gflops (gflops2=%lf,time=%lf)\n", 
          gflops2/time_2nd_vbatch, gflops2, time_2nd_vbatch );

   printf("overall %f gflops\n", (gflops1+gflops2)/(time_1st_vbatch + time_2nd_vbatch) );
#endif
     


 dmrg_free( BX_ );
}

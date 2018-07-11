#ifndef DMRG_VBATCH_H
#define DMRG_VBATCH_H 1

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

#ifndef ICEIL
#define ICEIL( x, n)  (( (x) + (n)-1 )/(n))
#endif

#ifndef indx2f
#define indx2f(i,j,lda)  (((i)-1) + ((j)-1)*(lda))
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern
void dmrg_init();



extern
void *dmrg_malloc( const size_t alloc_size );

extern
void  dmrg_free( void *a_ptr );

extern void
dmrg_dgetvector( const int n, 
                 const double *dx_src, const int incx, 
                       double *hy_dst, const int incy );

extern void
dmrg_dsetvector( const int n, 
                 const double *hx_src, const int incx, 
                       double *dy_dst, const int incy );

extern void
dmrg_dgetmatrix( const int m, const int n, 
                 const double *dA_src, const int ldda, 
                       double *hB_dst, const int ldb );

extern void
dmrg_dsetmatrix( const int m, const int n, 
                 const double *hA_src, const int lda, 
                       double *dB_dst, const int lddb );
        
extern
void dmrg_dgemm_vbatch( char transa_array[],
                        char transb_array[],
                        int  m_array[],
                        int  n_array[],
                        int  k_array[],
                        double alpha_array[],
                        double *a_array[], int    lda_array[],
                        double *b_array[], int    ldb_array[],
                        double beta_array[],
                        double *c_array[], int    ldc_array[],
                        int group_count, 
                        int group_size[] );

extern
void apply_Htarget_vbatch( 
                    int noperator,
                    int npatches, 
                    long left_patch_start_[],
                    long right_patch_start_[],
                    long xy_patch_start_[],
                    double Abatch_[], int ld_Abatch,
                    double Bbatch_[], int ld_Bbatch,
                    double X_[],
                    double Y_[]);
#ifdef __cplusplus
}
#endif

#endif

#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#include <getopt.h>
#include <string.h>


#ifdef _OPENMP
#include <omp.h>
#endif

#include "analysis.h"


#include "dmrg_vbatch.h"
#include "test_vbatch.h"

#ifdef USE_MAGMA
#include "cuda.h"
#include "magma_v2.h"
#endif

/*
 ---------------------------------------
 simple program to test vbatch
 ---------------------------------------
*/

int main(int argc, char *argv[])
{
 const int ialign = 32;
 const int idebug = 1;
 int noperator = 4;
 int left_size = 8;
 int max_keep_states = 4000;

#ifdef USE_MAGMA
  magma_init();
#endif
 dmrg_init();

 /*
  -----------------------------------------------------------
  ./test_vbatch -o noperator -n left_sites -m max_keep_states
  -----------------------------------------------------------
  */
 int opt = 0;
 while ((opt = getopt(argc,argv,"o:n:m:")) != -1) {
    switch(opt) {
     case 'n':
         left_size = atoi(optarg);
         break;
     case 'o':
         noperator = atoi(optarg);
         break;
     case 'm':
         max_keep_states = atoi(optarg);
         break;
     default: /* '?' */
         fprintf(stderr,"Usage: %s [-o noperator] [-n left_size] [-m max_states]\n",
            argv[0] );
         exit( EXIT_FAILURE );
    };
   };

  int right_size = left_size;
  int target_up = (left_size+right_size)/2;
  int target_down = target_up;
  int max_patches = (1 + target_up)*(1+target_down); 
/*
 ---------------------------
 assume left part is growing
 ---------------------------
*/
 const int keep_left_states = 4*max_keep_states;
 const int keep_right_states = max_keep_states;

 const int max_patches_dim = ialign * ICEIL( max_patches+1, ialign);

 int left_patch_size_[ max_patches_dim];
 int left_patch_up_[ max_patches_dim];
 int left_patch_down_[ max_patches_dim];

 int right_patch_size_[ max_patches_dim];
 int right_patch_up_[ max_patches_dim];
 int right_patch_down_[ max_patches_dim];

#define left_patch_size(i) left_patch_size_[(i)-1]
#define right_patch_size(i) right_patch_size_[(i)-1]

 int interaction_matrix_[ max_patches_dim * max_patches_dim ];
#define interaction_matrix(ipatch,jpatch) interaction_matrix_[indx2f(ipatch,jpatch,npatches)]

  {
   int nthreads = 1;
#ifdef _OPENMP
   #pragma omp parallel
   #pragma omp master
   { nthreads =  omp_get_num_threads(); }
#endif
   printf("using %d threads\n", nthreads );
   }


 int npatches = gen_patches_comb(
                     left_size,
                     right_size,
                     target_up,
                     target_down,
                     keep_left_states,
                     keep_right_states,

                     left_patch_size_,
                     right_patch_size_,

                     left_patch_up_,
                     left_patch_down_,

                     right_patch_up_,
                     right_patch_down_,

                     interaction_matrix_
                     );

/*
 ---------------------------
 estimate length of X vector
 ---------------------------
 */
 long xy_size = 0;
 {
 int ipatch = 0;
 for(ipatch=1; ipatch <= npatches; ipatch++) {
     int nrowX = right_patch_size(ipatch);
     int ncolX = left_patch_size(ipatch);
     xy_size += (nrowX * ncolX );
     };
 }
/*
 -------------------
 estimate total work
 -------------------
 */
 double total_flops_method1 = 0.0;
 double total_flops_method2 = 0.0;
 double total_flops = 0.0;
 {
 int ipatch = 0;
 int jpatch = 0;

 for(jpatch=1; jpatch <= npatches; jpatch++) {
 for(ipatch=1; ipatch <= npatches; ipatch++) {
    int has_work = interaction_matrix(ipatch,jpatch);
    if (has_work) {
      double flops_total = 0.0;
      double flops_method1 = 0.0;
      double flops_method2 = 0.0;
      int nrowA = left_patch_size(ipatch);
      int ncolA = left_patch_size(jpatch);
      int nrowB = right_patch_size(ipatch);
      int ncolB = right_patch_size(jpatch);
  
      cal_kron_flops( nrowA, nrowB, ncolA, ncolB, 
              &flops_total, &flops_method1,   &flops_method2);
  
      total_flops += flops_total;
      total_flops_method1 += flops_method1;
      total_flops_method2 += flops_method2;
      };
    };
    };
  };
  
 total_flops *= noperator;
 total_flops_method1 *= noperator;
 total_flops_method2 *= noperator;
    
    

 printf("left_size=%d, right_size=%d, target_up=%d, target_down=%d\n", 
         left_size,    right_size,    target_up,    target_down );
 printf("keep_left_states=%d, keep_right_states=%d, noperator=%d, xy_size=%ld\n",
         keep_left_states,    keep_right_states, noperator,xy_size );

 printf("npatches=%d\n", npatches );
 {
 double total_gflops = total_flops/(1000.0*1000.0*1000.0);
 double total_gflops_method1 = total_flops_method1 / (1000.0*1000.0*1000.0);
 double total_gflops_method2 = total_flops_method2 / (1000.0*1000.0*1000.0);

 printf("total_gflops=%lf, total_gflops_method1=%lf, total_gflops_method2=%lf\n",
         total_gflops,     total_gflops_method1,     total_gflops_method2 );

 };


 if (idebug >= 1) {
   int ipatch=0;
   for(ipatch=1; ipatch  <= npatches; ipatch++) {
     printf("ipatch=%d, left_patch_size=%d, right_patch_size=%d\n",
             ipatch, left_patch_size(ipatch), right_patch_size(ipatch) );
     };
  }

  


 double *Abatch_ = NULL;
 double *Bbatch_ = NULL;
 long *left_patch_start_ = NULL;
 long *right_patch_start_ = NULL;
 long *xy_patch_start_ = NULL;
 int ld_Abatch = 0;
 int ld_Bbatch = 0;
#define xy_patch_start(i) xy_patch_start_[(i)-1]
 
 setup_vbatch( noperator, npatches, 
               left_patch_size_, right_patch_size_,
               &left_patch_start_, &right_patch_start_, &xy_patch_start_, 
               &Abatch_, &ld_Abatch, &Bbatch_, &ld_Bbatch,
               interaction_matrix_
                );

#define Abatch(i,j)  Abatch_[indx2f(i,j,ld_Abatch)]
#define Bbatch(i,j)  Bbatch_[indx2f(i,j,ld_Bbatch)]


 int xy_size_dim = ialign * ICEIL( xy_size, ialign );
 double *X_  = (double *) dmrg_malloc( sizeof(double) * xy_size_dim );
 assert( X_ != NULL );

 double *Y_  = (double *) dmrg_malloc( sizeof(double) * xy_size_dim );
 assert( Y_ != NULL );

#define X(i) X_[ (i)-1]
#define Y(i) Y_[ (i)-1]

#define hX(i) hX_[ (i)-1]
#define hY(i) hY_[ (i)-1]
 
 {
#ifdef USE_SETVECTOR
  double hX_[xy_size_dim];
#else
  double *hX_ = X_;
#endif

  long i = 0;
  for(i=1; i <= xy_size; i++) {
     hX(i) = ((double) i )/( (double) xy_size );
     };

#ifdef USE_GETSET
   {
    const int incx = 1;
    const int incy = 1;
    dmrg_dsetvector( xy_size,
                   &(hX(1)), incx, 
                   X_, incy );
   };
#endif



  };

 
 {
 int itimes = 0;
 const int ntimes = 3;

 for(itimes=1; itimes <= ntimes; itimes++) {
   apply_Htarget_vbatch( noperator,npatches, 
                       left_patch_start_, right_patch_start_, xy_patch_start_,
                       Abatch_, ld_Abatch, Bbatch_,  ld_Bbatch, X_, Y_ );
    };
 }

#ifdef USE_GETSET
 double hY_[ xy_size_dim ];
 {
   const int incx = 1;
   const int incy = 1;
   dmrg_dgetvector( xy_size,
                    Y_,  incx,
                    &(hY(1)), incy );
 };
#else
 double *hY_ = Y_;
#endif

 /*
  * ---------------------------------
  * generate summary statistics for Y
  * ---------------------------------
  */
 {
  long i = 0;
  double Y_avg = 0;
  double Y_sd = 0;
  for(i=1; i <= xy_size; i++) {
    Y_avg += hY(i);
    };
  Y_avg = Y_avg/( (double) xy_size );
  for(i=1; i <= xy_size; i++) {
    Y_sd = Y_sd + (hY(i)-Y_avg)*(hY(i)-Y_avg);
    };
  Y_sd = sqrt( Y_sd );
  printf("Y_avg = %lf, Y_sd = %lf\n", Y_avg, Y_sd );
  };



#ifdef USE_MAGMA
  magma_finalize();
#endif

 exit(0);
 return(0);
}


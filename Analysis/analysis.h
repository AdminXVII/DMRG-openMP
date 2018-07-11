
#ifndef ANALYSIS_H

#define ANALYSIS_H 1

#include <math.h>

#ifndef MIN
#define MIN(x,y) (((x) < (y))?(x):(y))
#endif

#ifndef MAX
#define MAX(x,y) (((x) > (y))?(x):(y))
#endif

#ifndef MOD
#define MOD(x,y) ((x)%(y))
#endif

#ifndef ABS
#define ABS(x) (((x) > 0) ? (x) : (-(x)) )
#endif

#ifndef indx2f
#define indx2f(i,j,lda)  (((i)-1) + ((j)-1)*(lda))
#endif

#ifdef __cplusplus
extern "C" {
#endif

extern 
int gen_patches_comb( 
          int left_size,
          int right_size,
          int target_up, 
          int target_down, 
          int keep_left_states,
          int keep_right_states,

          int left_patch_size_[],
          int right_patch_size_[],

          int left_patch_up_[],
          int left_patch_down_[],
          int right_patch_up_[],
          int right_patch_down_[],

          int interaction_matrix_[]
          );

     
       
extern
void cal_kron_flops( 
           int nrowA, int nrowB, 
           int ncolA, int ncolB,
           double *pflops_total,
           double *pflops_method1,
           double *pflops_method2 );

#ifdef __cplusplus
}
#endif
 



   
#endif

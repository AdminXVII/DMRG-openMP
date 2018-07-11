#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "analysis.h"

/*
 ---------------------------------------
 simple program to test gen_patches_comb()
 ---------------------------------------
*/

int main()
{

 const int left_size = 8;
 const int right_size = 10;
 const int target_up = (left_size + right_size)/2;
 const int target_down = target_up;
 const int keep_left_states = 5000;
 const int keep_right_states = 4*keep_left_states;
 const int max_patches = (1 + target_up) * (1 + target_down);

 int left_patch_size_[ max_patches +1];
 int left_patch_up_[ max_patches +1];
 int left_patch_down_[ max_patches +1];

 int right_patch_size_[ max_patches +1];
 int right_patch_up_[ max_patches +1];
 int right_patch_down_[ max_patches +1];

#define left_patch_size(i) left_patch_size_[(i)-1]
#define right_patch_size(i) right_patch_size_[(i)-1]

 int interaction_matrix_[ (max_patches + 1)*(max_patches+1)];
#define interaction_matrix(ipatch,jpatch) interaction_matrix_[ indx2f(ipatch,jpatch,npatches) ]


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

 printf("left_size=%d, right_size=%d, target_up=%d, target_down=%d\n", 
         left_size,    right_size,    target_up,    target_down );
 printf("keep_left_states=%d, keep_right_states=%d\n",
         keep_left_states,    keep_right_states );

 printf("npatches=%d\n", npatches );
 {
 int ipatch=0;
 for(ipatch=1; ipatch  <= npatches; ipatch++) {
   printf("ipatch=%d, left_patch_size=%d, right_patch_size=%d\n",
           ipatch, left_patch_size(ipatch), right_patch_size(ipatch) );
   };
  }
 
 {
  int ipatch = 0;
  int jpatch = 0;
  int nnz = 0;

  for(jpatch=1; jpatch <= npatches; jpatch++) {
  for(ipatch=1; ipatch <= npatches; ipatch++) {
     int has_work = interaction_matrix(ipatch,jpatch);
     if (has_work) {
        nnz += 1;
        };
     };
     };
   printf("number of non-zeros in interaction_matrix is %d\n", nnz );
   
   }
  

 exit(0);
 return(0);
}


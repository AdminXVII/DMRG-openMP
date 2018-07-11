#ifndef TEST_VBATCH_H
#define TEST_VBATCH_H 1

#include <stdlib.h>
#include <math.h>
#include <assert.h>

extern
void setup_vbatch( int noperator,
                   int npatches,
                   int left_patch_size_[],
                   int right_patch_size_[],
                   long *left_patch_start_[],
                   long *right_patch_start_[],
                   long *xy_patch_start_[],
                   double *Abatch_[], int *ld_Abatch,
                   double *Bbatch_[], int *ld_Bbatch,
                   int interaction_matrix_[] );

#endif

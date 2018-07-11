/************************************************
 *
 * Created: 20th July, 2017
 * Updated: 27th July, 2017
 *
 ************************************************/

#include <errno.h>
#include <omp.h>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "Matrix.h"
#include "analysis.h"
#include "apply_Htarget.h"

//#define DEBUG_PATCHES 1

int main(int argc, char **argv) {
    char buf[256];
    int nrow, ncol, nnz;
    int ii, jj;
    double val;
    int n_c_rows = 0;
    int n_c_cols = 0;
    int vec_size = 0;
    auto C = new Block_Matrix_t;

    if (argc != 5) {
        printf(
            "Format: simulate_driver_generate SYSTEM_SIZE NUM_STATES "
            "SWEEP_LOCATION NUM_K_MAT\n");
        exit(1);
    }

    const int     is_grow_left = 1;
    const int max_keep_states = std::stoi(std::string(argv[2]));
    const int sys_size = std::stoi(std::string(argv[1]));
    int keep_left_states = 0, keep_right_states = 0;
    if (is_grow_left) {
        // left side is growing
         keep_left_states = 4*max_keep_states; 
         keep_right_states = max_keep_states; 
         }
     else {
        // right side is growing
        keep_left_states = max_keep_states; 
        keep_right_states = 4*max_keep_states;
        };
     
    const int left_size = std::stoi(std::string(argv[3]));
    const int NUM_K_MAT = std::stoi(std::string(argv[4]));
    const int right_size = sys_size - left_size;
    const int target_up = (left_size + right_size) / 2;
    const int target_down = target_up;
    const int max_patches = (1 + target_up) * (1 + target_down);

    int left_patch_size_[max_patches + 1];
    int left_patch_up_[max_patches + 1];
    int left_patch_down_[max_patches + 1];

    int right_patch_size_[max_patches + 1];
    int right_patch_up_[max_patches + 1];
    int right_patch_down_[max_patches + 1];

    int interaction_matrix_[(max_patches + 1) * (max_patches + 1)] = {-1};

#define left_patch_size(i) left_patch_size_[(i)-1]
#define right_patch_size(i) right_patch_size_[(i)-1]
#define interaction_matrix(ipatch,jpatch) interaction_matrix_[ indx2f((ipatch+1),(jpatch+1),npatches) ]


    int npatches = gen_patches_comb(
        left_size, right_size, target_up, target_down, keep_left_states,
        keep_right_states, left_patch_size_, right_patch_size_, left_patch_up_,
        left_patch_down_, right_patch_up_, right_patch_down_, &interaction_matrix_[0]);
/*    printf(
        "system_size = %d   left_states = %d   right_states = %d    left_size "
        "= %d   ",
        sys_size, keep_left_states, keep_right_states, left_size);
    printf(
        "right_size = %d   target_up = %d  target_down = %d   npatches = "
        "%d\n",
        right_size, target_up, target_down, npatches);
*/
    double lower_bound = 0.0;
    double upper_bound = 1.0;
    std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
    std::default_random_engine re;

    n_c_rows = npatches;
    n_c_cols = n_c_rows;
    C->cij.resize(n_c_rows);
    for (auto &a : C->cij) {
        a.resize(n_c_cols, nullptr);
    }

    std::vector<int> l_patch_size;
    l_patch_size.reserve(n_c_rows);
    std::vector<int> r_patch_size;
    r_patch_size.reserve(n_c_rows);

#ifdef DEBUG_PATCHES
    vec_size = 0;
    for (int ic = 0; ic < npatches; ic++) {
        vec_size += left_patch_size_[ic] * right_patch_size_[ic];
        for (int jc = 0; jc < npatches; jc++) {
            int nrowA, ncolA, nrowB, ncolB;
            nrowA = l_patch_size[ic] = left_patch_size_[ic];
            ncolA = l_patch_size[jc] = left_patch_size_[jc];
            nrowB = r_patch_size[ic] = right_patch_size_[ic];
            ncolB = r_patch_size[jc] = right_patch_size_[jc];

            printf("%d %d A %d %d B %d %d \n", ic+1, jc+1, nrowA, ncolA, nrowB,
                   ncolB);
        }
    }
    printf("vector_size = %d \n", vec_size);
    printf("Interaction matrix :\n");
    for (int ic = 0; ic < npatches; ic++){
        for (int jc = 0; jc < npatches; jc++){
            printf("%1d", interaction_matrix(ic, jc));
        }
        printf("\n");
    }
    exit(0);
#endif

    // for (auto line = 0; line < n_c_rows * n_c_cols; line++) {
    for (int ic = 0; ic < npatches; ic++) {
        vec_size += left_patch_size_[ic] * right_patch_size_[ic];
        for (int jc = 0; jc < npatches; jc++) {

            if (interaction_matrix(ic, jc) == 0){
                continue;
            }
            int nrowA, ncolA, nrowB, ncolB;
            nrowA = l_patch_size[ic] = left_patch_size_[ic];
            ncolA = l_patch_size[jc] = left_patch_size_[jc];
            nrowB = r_patch_size[ic] = right_patch_size_[ic];
            ncolB = r_patch_size[jc] = right_patch_size_[jc];

            //printf("%d %d A %d %d B %d %d \n", ic+1, jc+1, nrowA, ncolA, nrowB,
            //       ncolB);

            if (C->cij[ic][jc] == nullptr) {
                C->cij[ic][jc] = new CIJ_Elem_t;
            }
            CIJ_Elem celem = C->cij[ic][jc];
            celem->A.resize(NUM_K_MAT);
            celem->B.resize(NUM_K_MAT);
            for (int k = 0; k < NUM_K_MAT; k++) {
                auto newmat = new Matrix_t;
                newmat->nrow = nrowA;
                newmat->ncol = ncolA;
                newmat->is_dense = true;
                newmat->val.resize(nrowA * ncolA);
                for (int ii = 0; ii < nrowA * ncolA; ii++) {
                    newmat->val[ii] = ii * ii;
                }
                celem->A[k] = newmat;

                auto newmatB = new Matrix_t;
                newmatB->nrow = nrowB;
                newmatB->ncol = ncolB;
                newmatB->is_dense = true;
                newmatB->val.resize(nrowB * ncolB);
                for (int ii = 0; ii < nrowB * ncolB; ii++) {
                    newmatB->val[ii] = ii * ii;
                }
                celem->B[k] = newmatB;
            }
        }
    }
    std::vector<double> X, Y;
    X.reserve(vec_size);
    Y.reserve(vec_size);
    for (int i = 0; i < vec_size; i++) {
        X[i] = unif(re);
        Y[i] = 0.0;
    }

    double start, end;
    start = omp_get_wtime();
    int trials = 1;
    for (int num = 0; num < trials; num++) {
        //    CALLGRIND_START_INSTRUMENTATION;
        //    CALLGRIND_TOGGLE_COLLECT;
        apply_Htarget(*C, l_patch_size, r_patch_size, X, Y);
        //    CALLGRIND_TOGGLE_COLLECT;
        //    CALLGRIND_STOP_INSTRUMENTATION;
    }
    end = omp_get_wtime();
    std::cout << " Execution Time: " << ((double)(end - start) / trials)
              << std::endl;
}

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
#include <iostream>
#include <random>
#include <string>
#include <vector>
#include "Matrix.h"
#include "apply_Htarget.h"
#include "magma_kron.h"

int main(int argc, char **argv) {
  char buf[256];
  int nrow, ncol, nnz;
  int ii, jj;
  double val;
  int n_c_rows = 0;
  int n_c_cols = 0;
  int vec_size = 0;
  auto C = new Block_Matrix_t;
  int NUM_K_MAT = 0;
  double start, tall, tsetup, tcompute, tinit;

  if (argc != 3) {
    printf("Format: simulate_driver CIJ_FILE_NAME NUM_K_MAT\n");
    exit(1);
  }
  std::string cij_file_name(argv[1]);
  NUM_K_MAT = std::stoi(std::string(argv[2]));

  double lower_bound = 0.0;
  double upper_bound = 1.0;
  std::uniform_real_distribution<double> unif(lower_bound, upper_bound);
  std::default_random_engine re;

  std::ifstream ifs(cij_file_name, std::ifstream::in);
  ifs.getline(buf, 256);  // n_c_rows vec_size
  sscanf(buf, "%d %d", &n_c_rows, &vec_size);
  n_c_cols = n_c_rows;
  C->cij.resize(n_c_rows);
  for (auto &a : C->cij) {
    a.resize(n_c_cols, nullptr);
  }

  std::vector<int> l_patch_size;
  l_patch_size.reserve(n_c_rows);
  std::vector<int> r_patch_size;
  r_patch_size.reserve(n_c_rows);

  for (auto line = 0; line < n_c_rows * n_c_cols; line++) {
    int ic, jc, imatrix;
    char A, B;
    int nrowA, ncolA, nrowB, ncolB;
    ifs.getline(buf, 256);  // n_c_rows vec_size
    sscanf(buf, "%d %d %c %d %d %c %d %d", &ic, &jc, &A, &nrowA, &ncolA, &B,
           &nrowB, &ncolB);
    // printf("%d %d %c %d %d %c %d %d\n", ic, jc, A, nrowA, ncolA, B, nrowB,
    // ncolB);
    ic -= 1;
    jc -= 1;

    l_patch_size[ic] = nrowA;
    l_patch_size[jc] = ncolA;
    r_patch_size[ic] = nrowB;
    r_patch_size[jc] = ncolB;

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
        newmat->val[ii] = ii*ii;
      }
      celem->A[k] = newmat;

      auto newmatB = new Matrix_t;
      newmatB->nrow = nrowB;
      newmatB->ncol = ncolB;
      newmatB->is_dense = true;
      newmatB->val.resize(nrowB * ncolB);
      for (int ii = 0; ii < nrowB * ncolB; ii++) {
        newmatB->val[ii] = ii*ii;
      }
      celem->B[k] = newmatB;
    }
  }
  std::vector<double> X, Y;
  X.reserve(vec_size);
  Y.reserve(vec_size);
  for (int i = 0; i < vec_size; i++) {
    X[i] = unif(re);
    Y[i] = 0.0;
  }
  
  #if defined(MAGMA_CHECK_ERROR)
  std::vector < double > Ydup ( vec_size, 0.0);
  std::vector < double > work ( vec_size );
  #endif

  start = omp_get_wtime();
  tinit = start;
  magma_setup_kron_gpu(*C, l_patch_size, r_patch_size);
  tsetup = omp_get_wtime()-tinit;

  tinit = omp_get_wtime();
  int trials = 1;
  for (int num = 0; num < trials; num++) {
    magma_apply_Htarget_gpuv2(X, Y);
    #if defined(MAGMA_CHECK_ERROR)
    apply_Htarget(*C, l_patch_size, r_patch_size, X, Ydup);
    double c_neg_one =-1.0, kron_error=0.0, kron_error2=0.0, Ynorm;
    int ione = 1;
    daxpy_( &vec_size, &c_neg_one, &Ydup[0], &ione, &Y[0], &ione );
    Ynorm      = dlange_( "M", &vec_size, &ione, &Ydup[0], &vec_size, &work[0] );
    kron_error = dlange_( "M", &vec_size, &ione, &Y[0],    &vec_size, &work[0] );
    // compute elementwise | Y_magma - Y_default / Y_default |
    for (int j  = 0; j < vec_size; j++)
    {
        Y[j] = Y[j]/Ydup[j];
    }
    kron_error2 = dlange_( "M", &vec_size, &ione, &Y[0],    &vec_size, &work[0] );

    // |Y_magma - Y_default| / |Y_default|
    printf("===========================================\n");
    printf("    |Y_magma - Y_default| / |C_default|    \n");
    printf("    Magma_GPU-ref error is %10.5e\n",kron_error/Ynorm);
    printf("    | Y_magma - Y_default / Y_default |    \n");
    printf("    Magma_GPU-ref error is %10.5e\n",kron_error2);
    printf("===========================================\n");
    #endif
  }
  tcompute = omp_get_wtime()-tinit;

  magma_destroy_kron_gpu(*C, l_patch_size, r_patch_size);
  tall = omp_get_wtime()-start;

  int nthreads = 1;
  #if defined(MAGMA_CPU_USE_OMP)
  #pragma omp parallel
  #pragma omp master
      nthreads = omp_get_num_threads();
  #endif
 
  std::cout << " Magma Setup Time: " << ((double)(tsetup))
            << std::endl;
  std::cout << " Execution Time: " << ((double)(tcompute) / trials) <<  "   Gflops: " << ((flops_count*1e-9)/((double)(tcompute) / trials))  << " Number of threads: " << nthreads
            << std::endl;
  std::cout << " Total Time: " << ((double)(tall))
            << std::endl;

  
}

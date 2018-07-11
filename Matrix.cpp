#include "Matrix.h"


int Matrix_t::nnz() {
         return( (is_dense) ?   nrow*ncol : rowptr[nrow] );
     }

void Matrix_t::kron_mult( const char transA, 
                    const char transB, 
                    const Matrix_t A, 
                    const Matrix_t B, 
                    const double yin[], 
                    double xout[]) {
        
        const int nrow_A = A.nrow;
        const int ncol_A = A.ncol;
        const int nrow_B = B.nrow;
        const int ncol_B = B.ncol;

        const double *aval = &(A.val[0]);
        const double *bval = &(B.val[0]);
    
  
        const bool is_dense_A = A.is_dense;
        const bool is_dense_B = B.is_dense;

        const bool is_sparse_A = !is_dense_A;
        const bool is_sparse_B = !is_dense_B;
        

        if (is_dense_A) {
         if (is_dense_B) { 
           den_kron_mult(
                    transA, transB,
                    nrow_A, ncol_A, aval,
                    nrow_B, ncol_B, bval,
                    yin, xout );
           }
         else  {
           // B is sparse
           const int *browptr = &(B.rowptr[0]);
           const int *bcol = &(B.col[0]);

           den_csr_kron_mult(
                    transA, transB,
                    nrow_A, ncol_A, aval,
                    nrow_B, ncol_B, browptr, bcol, bval,
                    yin, xout );
              }
          }
        else {
           // A is sparse
           const int *arowptr = &(A.rowptr[0]);
           const int *acol = &(A.col[0]);

           if (is_dense_B) {
             csr_den_kron_mult(
                    transA, transB,
                    nrow_A, ncol_A, arowptr, acol, aval,
                    nrow_B, ncol_B, bval,
                    yin, xout );
             }
           else {
            // B is sparse
            const int *browptr = &(B.rowptr[0]);
            const int *bcol = &(B.col[0]);


            csr_kron_mult(
                    transA, transB,
                    nrow_A, ncol_A, arowptr, acol, aval,
                    nrow_B, ncol_B, browptr, bcol, bval,
                    yin, xout );
             };
          };
     } // kron_mult

          
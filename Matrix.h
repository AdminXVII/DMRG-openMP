#ifndef MATRIX_H
#define MATRIX_H

#include <vector>
#include "util.h"

class Matrix_t 
{
	public:
        std::vector < double > val;   // val[] holds the entries of matrix
        std::vector < int > rowptr;   // rowptr[i] holds the starting index in col[] where row "i" starts
        std::vector < int > col;      // column value of sparse matrix entry
        int nrow;
        int ncol;
        bool is_dense;                    // Is the matrix sparse / dense
        // note for dense matrix storage, then size(val) is nrow * ncol
        // consider using Fortran indexing   A(i,j) at  val[ (i) + (j)*nrow ]
        // to be compatible with BLAS
     int nnz(); 

     void kron_mult( const char transA, 
                    const char transB, 
                    const Matrix_t A, 
                    const Matrix_t B, 
                    const double yin[], 
                    double xout[]) ;
       
};

#endif

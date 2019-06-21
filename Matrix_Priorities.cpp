#include "Matrix.h"


int Matrix_t::nnz() {
         return( (is_dense) ?   nrow*ncol : rowptr[nrow] );
     }

#include "util.h"
int den_nnz( const int nrow_A, 
             const int ncol_A, 
             const double a_[])
{
#define A(ia,ja) a_[ (ia) + (ja)*nrow_A ]
/*
 * -------------------------
 * return number of nonzeros
 * matrix A in dense storage format
 * -------------------------
 */
  const double dzero = 0;
  int nnz_A = 0;

  int ja = 0;

  /*
   * -----------------------------
   * Note that  nnz_A is a reduction variable
   * -----------------------------
   */
  for(ja=0; ja < ncol_A; ja++) {
    int ia = 0;
    for(ia=0; ia < nrow_A; ia++) {
       int is_zero = (A(ia,ja) == dzero);
       nnz_A = (is_zero)? nnz_A : (nnz_A+1);
       };
     };

  return( nnz_A );
}
#undef A

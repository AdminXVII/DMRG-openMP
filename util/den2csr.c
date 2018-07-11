#include "util.h"
void den2csr( const int nrow_A, 
              const int ncol_A, 
              const double a_[], 

              const int max_nnz, 
                    int arowptr[], 
                    int acol[], 
                    double aval[] )
{
#define A(ia,ja) a_[(ia) + (ja)*nrow_A]

/*
 * -----------------------------
 * convert from dense matrix format to
 * compressed sparse row format
 * -----------------------------
 */

  const double dzero = 0;


  int nnz_row[nrow_A];
  
  int ia = 0;
  int ja = 0;


  for(ia=0; ia < nrow_A; ia++) {
    nnz_row[ia] = 0;
    };
  
  /*
   * -----------------------------------------------
   * first pass to count number of non zeros per row
   * -----------------------------------------------
   */

  for(ja=0; ja < ncol_A; ja++) {
  for(ia=0; ia < nrow_A; ia++) {
      double aij = A(ia,ja);
      int is_zero = (aij == dzero);
      if (!is_zero) {
           nnz_row[ia] += 1;
           };
      };
      };


  /*
   * --------------------------------------
   * check that there is sufficient storage
   * --------------------------------------
   */

   int  nnz = 0;
   for(ia=0; ia < nrow_A; ia++) {
     nnz += nnz_row[ia];
     };

   int isok = (nnz <= max_nnz);
   assert( isok );

 /*
  * --------------------------------
  * prefix sum to setup row pointers
  * --------------------------------
  */
  arowptr[0] = 0;
  for(ia=0; ia < nrow_A; ia++) {
    arowptr[ia+1] = arowptr[ia] + nnz_row[ia];
    };

 /*
  * -------------------------------------
  * second pass to fill in data structure
  * -------------------------------------
  */
 for(ia=0; ia < nrow_A; ia++) {
    nnz_row[ia] = 0;
    };


 for(ja=0; ja < ncol_A; ja++) {
 for(ia=0; ia < nrow_A; ia++) {
    double aij = A(ia,ja);
    int is_zero = (aij == dzero);
    if (!is_zero) {
       int ipos = arowptr[ia] + nnz_row[ia];
       acol[ipos] = ja;
       aval[ipos] = aij;
       nnz_row[ia] += 1;
       };
    };
    };

}
#undef A

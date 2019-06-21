#ifndef APPLY_HTARGET_H
#define APPLY_HTARGET_H

#include <vector>
#include "Matrix.h"

extern int Ntasks;

typedef Matrix_t *Matrix;

class CIJ_Elem_t
{
    public:
        std::vector < Matrix > A;     // List of A's - k items 
        std::vector < Matrix > B;     // List of B's - k items
};
typedef CIJ_Elem_t *CIJ_Elem;


class Block_Matrix_t
{
    public:
    	std::vector < std::vector< CIJ_Elem > > cij;   
};
typedef Block_Matrix_t *Block_Matrix;


typedef std::vector < int >  int_vec_t;
typedef std::vector < int_vec_t > int_vec_vec_t;

int apply_Htarget(Block_Matrix_t &CIJ, 
                  std::vector <int> &vsize,
                  std::vector <int> &vstart,
                  std::vector <double> &X,
                  std::vector <double> &Y);

void kron_mult( const char transA, const char transB,
				const Matrix A, const Matrix B,
				const double yin[], double xout[] );

#endif


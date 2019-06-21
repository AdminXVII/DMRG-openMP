/************************************************
 * Author: Arghya Chatterjee
 * Created: 26th Jan, 2017
 * Updated: 27th July, 2017
 * Multi-Level Parallel loops on-and-off
 ************************************************/

#include "apply_Htarget.h"
#include <iostream>
#include <array>
#include <vector>
#include <omp.h>
#include <algorithm>

#define tH 1000
#define NITS 10 
#define NThreads 40
#define PRIOTH 7000

void kron_mult( const char transA, const char transB,
				const Matrix A, const Matrix B,
				const double yin[], double xout[] ){
	const int nrow_A = A->nrow;
    const int ncol_A = A->ncol;
    const int nrow_B = B->nrow;
    const int ncol_B = B->ncol;

    const double *aval = &(A->val[0]);
    const double *bval = &(B->val[0]);
    
    const bool is_dense_A = A->is_dense;
    const bool is_dense_B = B->is_dense;

    if (is_dense_A) {
    	if (is_dense_B) { 
        	den_kron_mult(transA, transB,
                    nrow_A, ncol_A, aval,
                    nrow_B, ncol_B, bval,
                    yin, xout );
        }
		else {
        	// B is sparse
            const int *browptr = &(B->rowptr[0]);
            const int *bcol = &(B->col[0]);
            den_csr_kron_mult(transA, transB,
                    nrow_A, ncol_A, aval,
                    nrow_B, ncol_B, browptr, bcol, bval,
                    yin, xout );
        }
	}
    else {
    	// A is sparse
        const int *arowptr = &(A->rowptr[0]);
        const int *acol = &(A->col[0]);

        if (is_dense_B) {
        	csr_den_kron_mult(transA, transB,
                    nrow_A, ncol_A, arowptr, acol, aval,
                    nrow_B, ncol_B, bval,
                    yin, xout );
        }
        else {
        	// B is sparse
            const int *browptr = &(B->rowptr[0]);
            const int *bcol = &(B->col[0]);
            csr_kron_mult(transA, transB,
                    nrow_A, ncol_A, arowptr, acol, aval,
                    nrow_B, ncol_B, browptr, bcol, bval,
                    yin, xout );
        }
	}
}


int apply_Htarget(Block_Matrix_t &CIJ, std::vector<int> &vsize,
				  std::vector<int> &vstart, std::vector < double > &X,
                  std::vector < double > &Y){
    int npatches = CIJ.cij[0].size();

	double* Y_ptr = &Y[0];
    double* X_ptr = &X[0];
    char* sentinel = new char[npatches]();
    #pragma omp parallel
    #pragma omp single
    { // start parallel region for iPatch
       for(int its = 0; its < NITS; its++){ 
        for(int ipatch = 0; ipatch < npatches; ipatch++)
        {
            int i1 = vstart[ipatch];
            int i2 = i1 + vsize[ipatch];
            int fine_grain = (vsize[ipatch] <= tH);
            int coarse_grain = (vsize[ipatch] >= tH);
            int prio_rosa = 10*coarse_grain + 2; //O 12 o 2
            #pragma omp task firstprivate(ipatch, i1, i2, fine_grain) depend(inout: sentinel[ipatch]) final(fine_grain) priority(prio_rosa) 
            {
            for(int jpatch = 0; jpatch < npatches; jpatch++)
            {
            	int j1 = vstart[jpatch];
            	int j2 = j1 + vsize[jpatch];
                int size_list_k = CIJ.cij[ipatch][jpatch] == nullptr ? 0 :  
                CIJ.cij[ipatch][jpatch] -> A.size();
                
                for(int k = 0; k < size_list_k; k++){
                    Matrix Ak = CIJ.cij[ipatch][jpatch] -> A[k];
                    Matrix Bk = CIJ.cij[ipatch][jpatch] -> B[k];
                    if ( Ak -> nnz() && Bk -> nnz() ){ //Lazy evaluation!
                        int DIAGONAL = (ipatch == jpatch); //Better not use branches
                        if(fine_grain){
                        	//#pragma omp task depend(inout: Y_ptr[i1:i2]) depend(in: X_ptr[j1:j2]) firstprivate(Ak, Bk) priority(0)
                        		{
                        		kron_mult('n','n', Ak, Bk, &X_ptr[j1], &Y_ptr[i1]); 
                        		}
                        }
                        else{
                        	double** buffer = new double*;
							
                        	#pragma omp task depend(out: buffer) depend(in: X_ptr[j1:j2]) firstprivate(Ak, Bk, buffer, ipatch) shared(X_ptr,vsize) priority(0) 
                        	{
                        		double* Y_return = new double[vsize[ipatch]]();
                        		buffer[0] = Y_return;
                        		kron_mult('n','n', Ak, Bk, &X_ptr[j1], Y_return);
                        	}
                        	#pragma omp task depend(inout: Y_ptr[i1:i2]) depend(in: buffer) firstprivate (i1, i2, buffer) shared(Y_ptr) priority(1)
                        	{
                            	double* Y_return = buffer[0];
                            	int ilocal = 0;
                            	int i;
                                #pragma unroll 8
                                for(i=i1; i < i2; i++)
                                	Y_ptr[i] += Y_return[ilocal++];
                            	delete[] Y_return;
                            	delete[] buffer;
                        	}
						}	
                    }//close has work 
                } // end of k loop
                //#pragma omp taskwait
            } // end of jpatch loop
            #pragma omp taskwait
            } // end of task
        } // end of ipatch
        }
    } //close parallel single 
    delete[] sentinel;
    std::cout<<"Done ApplyHTarget"<<"\n\n"; 
    return 1;
}// end apply_Htarget

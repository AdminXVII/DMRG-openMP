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
#define NBUFF 100 
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


int apply_Htarget(Block_Matrix_t &CIJ, int_vec_t &Lindex_Patch,
                  int_vec_t &Rindex_Patch, std::vector < double > &X,
                  std::vector < double > &Y){
    
    int npatches = CIJ.cij[0].size();
    std::vector < int > vsize(npatches);
    std::vector < int > vstart(npatches);
    
    int nrowL;       // No of rows for the Lindex
    int nrowR;       // No of rows for the Rindex
    int ip = 0;      // ip: Points to start of each patch
    int DIAGONAL = 0;
     
    for(int ipatch = 0; ipatch < npatches; ipatch++){
        nrowL = Lindex_Patch[ipatch];
        nrowR  = Rindex_Patch[ipatch];
        vsize[ipatch] = nrowL * nrowR;
	}
    
    for(int ipatch = 0; ipatch < npatches; ipatch++){
        vstart[ipatch] = ip;
        ip = ip + vsize[ipatch];
    }
    
    double* buffers[NBUFF];
    int next = 0;

	double* Y_ptr = &Y[0];
    double* X_ptr = &X[0];

    #pragma omp parallel 
    #pragma omp single
    { // start parallel region for iPatch
        
        for(int ipatch = 0; ipatch < npatches; ipatch++)
        {
            int i1 = vstart[ipatch];
            int i2 = i1 + vsize[ipatch];

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
                        DIAGONAL = (ipatch == jpatch); //Better not use branches
                        if(vsize[ipatch] <= tH){
                        	#pragma omp task depend(inout: Y_ptr[i1:i2]) depend(in: X_ptr[j1:j2]) firstprivate(Ak, Bk)  priority(0)
                        		{
                        		//Matrix MA = (CIJ.cij[ipatch][jpatch] -> A[k]);
                        		//Matrix MB = (CIJ.cij[ipatch][jpatch] -> B[k]);
                        		kron_mult('n','n', Ak, Bk, &X_ptr[j1], &Y_ptr[i1]); 
                        		}
                        }
                        else{
                        	int mybuff;
                        	mybuff = next = (next+1)%NBUFF;
							int prio = (vsize[ipatch] > PRIOTH);
                        	#pragma omp task depend(inout: buffers[mybuff]) depend(in: X_ptr[j1:j2]) firstprivate(Ak, Bk, mybuff, ipatch) shared(X_ptr, vsize) priority(prio)
                        	{
                        		double* Y_return = new double[vsize[ipatch]]();
                        		buffers[mybuff] = Y_return;
                        		//Matrix MA = (CIJ.cij[ipatch][jpatch] -> A[k]);
                        		//Matrix MB = (CIJ.cij[ipatch][jpatch] -> B[k]);
                        		kron_mult('n','n', Ak, Bk, &X_ptr[j1], Y_return); 
                        	}
                        	#pragma omp task depend(inout: Y_ptr[i1:i2], buffers[mybuff]) firstprivate (i1, i2, mybuff) shared(Y_ptr) priority(10)
                        	{
                            	double* Y_return = buffers[mybuff];
                            	int ilocal = 0;
                            	//int lim = i2 - 7;
                            	int i;
                            	/*for(i = i1; i < lim; i+=8){
                                	Y_ptr[i] += Y_return[ilocal];
                                	Y_ptr[i+1] += Y_return[ilocal+1];
                                	Y_ptr[i+2] += Y_return[ilocal+2];
                                	Y_ptr[i+3] += Y_return[ilocal+3];
                                	Y_ptr[i+4] += Y_return[ilocal+4];
                                	Y_ptr[i+5] += Y_return[ilocal+5];
                                	Y_ptr[i+6] += Y_return[ilocal+6];
                                	Y_ptr[i+7] += Y_return[ilocal+7];
                                	ilocal += 8;
                                }*/
                                #pragma unroll 8
                                for(i=i1; i < i2; i++)
                                	Y_ptr[i] += Y_return[ilocal++];
                            	delete[] Y_return;
                            	buffers[mybuff]=NULL;
                        	}
						}	
                    }//close has work 
                } // end of k loop
            } // end of jpatch loop
        } // end of ipatch
    } //close parallel single 
    std::cout<<"Done ApplyHTarget"<<"\n\n"; 
    return 1;
}// end apply_Htarget

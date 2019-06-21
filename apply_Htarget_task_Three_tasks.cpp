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

#define tH 300

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

    //omp_set_num_threads( 40 );
    
    double* buffers[50];
    int next = 0;
    double* Y_ptr = &Y[0];
    double* X_ptr = &X[0];

    #pragma omp parallel 
    #pragma omp single
    { // start parallel region for iPatch
        
        for(int ipatch = 0; ipatch < npatches; ipatch++)
        {
            int i1, i2;
            i1 = vstart[ipatch];
            i2 = i1 + vsize [ipatch]; 
            
            int jpatch_start = 0;
            int jpatch_end = npatches;
            
            for(int jpatch = jpatch_start; jpatch < jpatch_end; jpatch++)
            {
                int j1, j2;
                j1 = vstart[jpatch];
                j2 = j1 + vsize[jpatch]; 
                int size_list_k = CIJ.cij[ipatch][jpatch] == nullptr ? 0 :  
                CIJ.cij[ipatch][jpatch] -> A.size();
                
                for(int k = 0; k < size_list_k; k++){
                    Matrix Ak = CIJ.cij[ipatch][jpatch] -> A[k];
                    Matrix Bk = CIJ.cij[ipatch][jpatch] -> B[k];
                    if ( Ak -> nnz() && Bk -> nnz() ){ //Lazy evaluation!
                        DIAGONAL = (ipatch == jpatch); //Better not use branches
                        if(vsize[ipatch] <= tH){
                            #pragma omp task depend(inout: Y_ptr[i1:i2]) depend(in: X_ptr[j1:j2], Ak, Bk)
                            {
                                Ak->kron_mult ('n','n', *Ak, *Bk, &X_ptr[j1], &Y_ptr[i1]); 
                            }
                        }
                        else{
                            int mybuff;
                            mybuff = next = (next+1)%50;
                            double* Y_return = new double[vsize[ipatch]]();
                            buffers[next] = &Y_return[0];
                            #pragma omp task depend(inout: Y_return[0:vsize[ipatch]]) depend(in: X_ptr[j1:j2], Ak, Bk)
                            Ak->kron_mult ('n','n', *Ak, *Bk, &X_ptr[j1], &Y_return[0]); 
                            int ilocal=0;
                            #pragma omp task depend(inout: Y_ptr[i1:i2]) depend(in: Y_return[0:vsize[ipatch]]) firstprivate (ilocal, i1, i2)
                            {
                                for(int i = i1; i < i2; i++)
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

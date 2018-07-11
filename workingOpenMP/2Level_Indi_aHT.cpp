/************************************************
 * Author: Arghya Chatterjee
 * Created: 26th Jan, 2017
 * Updated: 17th May, 2017
 *
 ************************************************/

#include "apply_Htarget.h"
#include <iostream>
#include <array>
#include <vector>
#include <omp.h>
#include <algorithm>


int apply_Htarget(Block_Matrix_t &CIJ, int_vec_t &Lindex_Patch,
                  int_vec_t &Rindex_Patch, std::vector < double > &X,
                  std::vector < double > &Y){
    int npatches = CIJ.cij[0].size();

    std::vector < int > vsize(npatches);
    std::vector < int > vstart(npatches);
 
    /**
      *  Calculating the size of the patches 
      *  and row and column indeces for X and Y
      *  matrices.
    **/
    
    int nrowL;       // No of rows for the Lindex
    int nrowR;       // No of rows for the Rindex
    int ip = 0;      // ip: Points to start of each patch
    int diagonal = 0;

    /**
      *  Calculate length of the X vector to be used for   
      *  each corresponding patch
      *  Version1 : Calculated from LIndex and RIndex
    **/

      
    for(int ipatch = 0; ipatch < npatches; ipatch++){
        nrowL = Lindex_Patch[ipatch];
        nrowR  = Rindex_Patch[ipatch];
        vsize[ipatch] = nrowL * nrowR;
    }
    

    /**
      *  Calculate length of the X vector to be used for   
      *  each corresponding patch
      *  Version2 : Calculated using patch size for each CIJ patch
    **/

    /** 
      *  for(int ipatch = 0; ipatch < npatches; ipatch++){
      *      nrowR = 0;
      *      nrowL = 0;
      *      for(int jpatch = 0; jpatch < npatches; jpatch++){
      *          if((CIJ.cij[ipatch][jpatch] -> A).size() >= 1){
      *              nrowL = CIJ.cij[ipatch][jpatch] -> A[0] -> ncol;
      *              nrowR = CIJ.cij[ipatch][jpatch] -> B[0] -> ncol;
      *          }     
      *      }
      *      vsize[ipatch] = nrowL * nrowR;
      *  }
    **/

    for(int ipatch = 0; ipatch < npatches; ipatch++){
        vstart[ipatch] = ip;
        ip = ip + vsize[ipatch];
    }

    #pragma omp parallel for
    for(int ipatch = 0; ipatch < npatches; ipatch++){
        int i1, i2;
        i1 = vstart[ipatch];
        i2 = i1 + vsize [ipatch]; 
        
        // needs to be private for each iteration 
        std::vector < double > YI ( vsize[ipatch] , 0.0);
        
        /**
          * Symmetry : Access only upper triangular matrix
        **/
        #ifdef use_symmetry
               int jpatch_start = 0;
               int jpatch_end = ipatch;
        #else
               int jpatch_start = 0;
               int jpatch_end = npatches;
        #endif

        #ifdef use_symmetry
               int size_XI = i2 - i1;
               std::vector < double > XI( size_XI );
               for(int i = i1; i < i2; i++){
               XI[i] = X[i];
               }
        #endif

        #pragma omp declare reduction(vec_double_add : std::vector < double >:\
                       std::transform(omp_in.begin(), omp_in.end(), \ 
                                      omp_out.begin(), omp_out.begin(), \
                                      std::plus< double > ())) \
                          initializer(omp_priv (omp_orig))
                          
        #pragma omp parallel reduction(vec_double_add : YI)
        {               
            #pragma omp for          
            for(int jpatch = jpatch_start; jpatch < jpatch_end; jpatch++){
                int j1, j2;
                j1 = vstart[jpatch];
                j2 = j1 + vsize[jpatch]; 

                int size_XJ = j2 - j1;
                std::vector < double > XJ (size_XJ);
                for(int j = j1; j < j2; j++){
                    XJ[j - j1] = X[j];
                }

                std::vector < double > YIJ ( vsize[ipatch], 0.0 );

                int size_list_k = CIJ.cij[ipatch][jpatch] == nullptr ? 0 : 
                                  CIJ.cij[ipatch][jpatch] -> A.size();

                for(int k = 0; k < size_list_k; k++){
                    std::vector < double > Y_return ( vsize[ipatch], 0.0 );
                    Matrix Ak = CIJ.cij[ipatch][jpatch] -> A[k];
                    Matrix Bk = CIJ.cij[ipatch][jpatch] -> B[k];
                    
                    int has_work = ( Ak -> nnz() && Bk -> nnz() );
                    if(!has_work){
                       continue;
                    }
                    if(ipatch == jpatch){
                       diagonal = 1;  
                    }
                    fflush(stdout);
                    #if defined (use_symmetry) && (diagonal)
                        Ak->kron_mult('n','n', *Ak, *Bk, &XJ[0], &Y_return[0]); 
                        for(int i = 0; i < vsize[ipatch]; i++){
                            YIJ[i] += Y_return[i];
                        }
                        Ak->kron_mult('t','t', *Ak, *Bk, &XI[0], &Y_return[0]);
                        for(int j = j1; j < j2; j++){
                            Y[j] += Y_return [j - j1];
                        }
                    #else
                        Ak->kron_mult ('n','n', *Ak, *Bk, &XJ[0], &Y_return[0]); 
                        for(int i = 0; i < vsize[ipatch]; i++){
                               YIJ[i] += Y_return[i];
                           }
                    #endif
                }// end of k block
                for(int i = 0; i< vsize[ipatch]; i++){
                    YI[i] += YIJ[i];            
                }
            }// end of jpatch        
        } // end parallel region
        for(int i = i1; i < i2; i++){
            Y[i] = YI[i - i1];
        }
    }     // end of ipatch
    return 1;

}     // end apply_Htarget



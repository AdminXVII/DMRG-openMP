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


#define NUM_THREADS_L0 
#define NUM_THREADS_L1 
#define NUM_THREADS_L2 

#define BIND_L0 spread
#define BIND_L1 close
#define BIND_L2 close

//#define PARALLEL_IPATCH

//choose mode:: (to use either mode, PARALLEL_IPATCH must be turned on)
//#define PARALLEL_IPATCH_WORKSHARING
//#define PARALLEL_IPATCH_TASKING

//#define PARALLEL_JPATCH 

//choose mode:: (to use either mode, PARALLEL_JPATCH must be turned on)
//#define PARALLEL_JPATCH_WORKSHARING
//#define PARALLEL_JPATCH_TASKING

//#define PARALLEL_KPATCH 
//#define CALCULATE_WORK

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
    int DIAGONAL = 0;

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


    int TThreads = omp_get_max_threads();

    struct Element
     {
        int Thread_L0;
        int Thread_L1;
        int CPU_ID;
        int work;
     }; 

    struct Element TUsage[npatches][npatches];

#ifdef PARALLEL_IPATCH
   #pragma omp parallel num_threads (NUM_THREADS_L0) proc_bind(BIND_L0)
#endif
    { // start parallel region for iPatch
 
#ifdef PARALLEL_IPATCH_WORKSHARING
             #pragma omp for schedule (dynamic,1)
#endif    
       
             for(int ipatch = 0; ipatch < npatches; ipatch++)
             {

#ifdef PARALLEL_IPATCH_TASKING
                    #pragma omp single
                    #pragma omp task 

#endif 
                    {  //open I tasking

                    int CThread_L_Zero = omp_get_thread_num(); 
                    int i1, i2;
                    i1 = vstart[ipatch];
                    i2 = i1 + vsize [ipatch]; 
                    
                    std::vector < double > YI ( vsize[ipatch] , 0.0);
                    
                    /**
                      * Symmetry : Access only upper triangular matrix
                    **/
#ifdef USE_SYMMETRY
                           int jpatch_start = 0;
                           int jpatch_end = ipatch;
#else
                           int jpatch_start = 0;
                           int jpatch_end = npatches;
#endif


#ifdef USE_SYMMETRY
                           int size_XI = i2 - i1;
                           std::vector < double > XI( size_XI );
                           for(int i = i1; i < i2; i++){
                           XI[i - i1] = X[i];
                           }
#endif

#ifdef PARALLEL_JPATCH                                                        
                    #pragma omp declare reduction(vec_double_add : std::vector < double >:\
                                   std::transform(omp_in.begin(), omp_in.end(),\
                                                  omp_out.begin(), omp_out.begin(),\
                                                  std::plus< double > ())) \
                                      initializer(omp_priv (omp_orig))
                    #pragma omp parallel num_threads(NUM_THREADS_L1)\
                                           proc_bind(BIND_L1)\
                                           reduction(vec_double_add : YI)
                  
#endif
                    {   // starting parallel region for jPatch            

#ifdef PARALLEL_JPATCH_WORKSHARING                                                        
#pragma omp for schedule (dynamic,1)
#endif                      
                        for(int jpatch = jpatch_start; jpatch < jpatch_end; jpatch++)
                        {
                            int CThread_L_One = omp_get_thread_num();

#ifdef PARALLEL_JPATCH_TASKING
                            #pragma omp single
                            #pragma omp task 
#endif 
                            {  //open J tasking


#ifdef CALCULATE_WORK                          
                            TUsage[ipatch][jpatch].CPU_ID = sched_getcpu();
                            TUsage[ipatch][jpatch].Thread_L0 = CThread_L_Zero;
                            TUsage[ipatch][jpatch].Thread_L1 = CThread_L_One;                
#endif

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

#ifdef CALCULATE_WORK
                            TUsage[ipatch][jpatch].work = size_list_k;
#endif

#ifdef PARALLEL_KPATCH                                                        
                            #pragma omp parallel num_threads(NUM_THREADS_L2)\
                                                    proc_bind(BIND_L2)\
                                                    reduction(vec_double_add: YIJ)
#endif
                            { // starting parallel region for kPatch

#ifdef PARALLEL_KPATCH                                                        
                                #pragma omp for schedule(dynamic,1)
#endif
                                for(int k = 0; k < size_list_k; k++)
                                {
                                    std::vector < double > Y_return ( vsize[ipatch], 0.0 );
                                    Matrix Ak = CIJ.cij[ipatch][jpatch] -> A[k];
                                    Matrix Bk = CIJ.cij[ipatch][jpatch] -> B[k];
                                    
                                    int has_work = ( Ak -> nnz() && Bk -> nnz() );
                                    if(!has_work){
                                       continue;
                                    }
                                    if(ipatch == jpatch){
                                       DIAGONAL = 1;  
                                    }
#if defined (USE_SYMMETRY) && (DIAGONAL)
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
                                } // end of k block
                            } // end of parallel region K


                            } //close J tasking 

                            for(int i = 0; i< vsize[ipatch]; i++){
                                YI[i] += YIJ[i];            
                            }

                       } //end of jpatch

                    } // end of parallel region J


                    for(int i = i1; i < i2; i++){
                        Y[i] = YI[i - i1];
                    }

                } //close I tasking

            } // end of ipatch

    } //end of parallel region I

#ifdef CALCULATE_WORK
     for (int i = 0; i < npatches; i++){
         for(int j = 0; j < npatches; j++){
            std::cout<< std::flush;

             std::cout<< "Patch: [ "<< i << " , " << j << " ] ::  " 
                     // << TUsage[i][j].Thread_L0 << "  ==>  "
                     // << TUsage[i][j].Thread_L1 << "  ==>  "
                     // << TUsage[i][j].CPU_ID << std::endl;
                    << TUsage[i][j].work << std::endl;
         }
         std::cout<<"\n\n";
    }
#endif

  std::cout<<"Done ApplyHTarget"<<"\n\n"; 
 
    return 1;

}
     // end apply_Htarget



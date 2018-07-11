/*******************************************************
 * Author: Arghya Chatterjee & Jun Shirako (HabaneroC++)
 * Created: 26th Jan, 2017
 * Updated: 28th August, 2017
 *
 *******************************************************/

#include "apply_Htarget.h"
#include <iostream>
#include <array>
#include <vector>
#include <omp.h>
#include <algorithm>

#include "accumulator.h"

// Undefine this to see ideal performance without accumulators (but incorrect output)
#define USE_ACCUMULATOR

/*
#define NUM_THREADS_L0 10
#define NUM_THREADS_L1 8
#define NUM_THREADS_L2 2
#define BIND_L0 spread
#define BIND_L1 close
#define BIND_L2 close
*/

#define NUM_THREADS_L0 24
#define BIND_L0 spread
#define BIND_L1 close
#define BIND_L2 close 

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


    int TThreads = omp_get_max_threads();

    struct Element
     {
        int Thread_L0;
        int Thread_L1;
        int CPU_ID;
     }; 

    struct Element TUsage[npatches][npatches];

#ifdef USE_ACCUMULATOR
    std::vector<accumulator *> accs(npatches);
    for (int i = 0; i < npatches; i++)
        accs[i] = new accumulator();
#endif

    #pragma omp parallel num_threads (NUM_THREADS_L0) proc_bind(BIND_L0)
    // Note: also set IBM smprt options, e.g., when using 24 workers (cores):
    //  $> export XLSMPOPTS=PARTHDS=24:SPINS=0:YIELDS=0:STARTPROC=0:STRIDE=8
    {
  //         #pragma omp for
           #pragma omp single
           {
             for(int ipatch = 0; ipatch < npatches; ipatch++){
                #pragma omp task shared (TUsage)
                {
                    int CThread_L_Zero = omp_get_thread_num(); 
                    int i1, i2;
                    i1 = vstart[ipatch];
                    i2 = i1 + vsize [ipatch]; 
                    
                    // needs to be private for each iteration 
                    // std::vector < double > YI ( vsize[ipatch] , 0.0);
		    // Directly added to accumulator
                    
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
                           XI[i - i1] = X[i];
                           }
                    #endif

                  /*
                    #pragma omp declare reduction(vec_double_add : std::vector < double >:\
                                   std::transform(omp_in.begin(), omp_in.end(),\
                                                  omp_out.begin(), omp_out.begin(),\
                                                  std::plus< double > ())) \
                                      initializer(omp_priv (omp_orig))
                                      
                   #pragma omp parallel num_threads(NUM_THREADS_L1)\
                                           proc_bind(BIND_L1)\
                                           reduction(vec_double_add : YI)
                  */

                    //#pragma omp taskgroup                   
                    {               

                      //define a reduction variable (i instances of variables)

                        //#pragma omp for 
                        for(int jpatch = jpatch_start; jpatch < jpatch_end; jpatch++){
                            #pragma omp task shared (TUsage)//, reduction variable)
                            {

                            int CThread_L_One = 0; //omp_get_thread_num();
                            TUsage[ipatch][jpatch].CPU_ID = sched_getcpu();
                            TUsage[ipatch][jpatch].Thread_L0 = CThread_L_Zero;
                            TUsage[ipatch][jpatch].Thread_L1 = CThread_L_One;                

                            int j1, j2;
                            j1 = vstart[jpatch];
                            j2 = j1 + vsize[jpatch]; 

                            int size_XJ = j2 - j1;
                            std::vector < double > XJ (size_XJ);
                            for(int j = j1; j < j2; j++){
                                XJ[j - j1] = X[j];
                            }

                            // std::vector < double > YIJ ( vsize[ipatch], 0.0 );
			    // Directly added to accumulator

                            int size_list_k = CIJ.cij[ipatch][jpatch] == nullptr ? 0 :  
                                              CIJ.cij[ipatch][jpatch] -> A.size();

                            // #pragma omp parallel num_threads(NUM_THREADS_L2)\
                            //                        proc_bind(BIND_L2)\
                            //                        reduction(vec_double_add: YIJ)

                            {
                                //#pragma omp for
                                for(int k = 0; k < size_list_k; k++){
                                   // #pragma omp task 
                                    {
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
                                    #if defined (use_symmetry) && (diagonal)
                                        Ak->kron_mult('n','n', *Ak, *Bk, &XJ[0], &Y_return[0]); 
                                        // for(int i = 0; i < vsize[ipatch]; i++){
                                        //     YIJ[i] += Y_return[i];
                                        // }
#ifdef USE_ACCUMULATOR
					accs[ipatch]->put(Y_return);
#endif
                                        Ak->kron_mult('t','t', *Ak, *Bk, &XI[0], &Y_return[0]);
                                        for(int j = j1; j < j2; j++){
                                            Y[j] += Y_return [j - j1];
                                        }
                                    #else
                                        Ak->kron_mult ('n','n', *Ak, *Bk, &XJ[0], &Y_return[0]); 
                                        // for(int i = 0; i < vsize[ipatch]; i++){
                                        //     YIJ[i] += Y_return[i];
                                        // }
#ifdef USE_ACCUMULATOR
					accs[ipatch]->put(Y_return);
#endif
                                    #endif
                                    }//end of task -- l2
                                }// end of k block
                            } // end parallel region -- K

                            // for(int i = 0; i< vsize[ipatch]; i++){
                            //     YI[i] += YIJ[i];            
                            // }

                          }//inner task -- l1

                        }// end of jpatch   

                    } // task group -- J

                    // for(int i = i1; i < i2; i++){
                    //     Y[i] = YI[i - i1];
                    // }
                    // To be updated via accumulator after parallel region

               }//end of outer task -- l0


            }     // end of ipatch

       }//outermost single     

    }//end parallel region -- I

#ifdef USE_ACCUMULATOR
    for (int ipatch = 0; ipatch < npatches; ipatch++) {
        accs[ipatch]->finish();

        int i1 = vstart[ipatch];
        int i2 = i1 + vsize[ipatch];
        std::vector<double> YI = accs[ipatch]->get();
        if (!YI.empty()) {
            for(int i = i1; i < i2; i++)
            Y[i] = YI[i - i1];
        }
    }
#endif
/*
     for (int i = 0; i < npatches; i++){
         for(int j = 0; j < npatches; j++){
            std::cout<< std::flush;

             std::cout<< "Patch: [ "<< i << " , " << j << " ] ::  " 
                      << TUsage[i][j].Thread_L0 << "  ==>  "
                      << TUsage[i][j].Thread_L1 << "  ==>  "
                      << TUsage[i][j].CPU_ID << std::endl;
         }
         std::cout<<"\n\n";
    }

  std::cout<<"Done ApplyHTarget"<<"\n\n"; 
*/
 
    return 1;

}     // end apply_Htarget



================================================================================
Project Title :
================================================================================
Co-Designing DMRG++ miniapp for Exascale Machines


================================================================================
List of files :
================================================================================
seq_apply_Htarget.cpp   -- Y = CIJ.(Ak kron Bk).X 
pl_reduce_2Level_aHT    -- Parallel Version (2 Levels - iPatch and jPatch)
apply_Htarget.cpp       -- Will be used as a running update for OpenMP
3Level_Indi_aHT.cpp     -- 3 Levels OpenMP (Individual layers of parallelism)
simulate_driver.cpp     -- Traverse through FakeCIJ shape for performance analysis
traverse.cpp            -- Traverse through CIJ dump and execute apply_Htarget
affinity.cpp            -- Module file for OmpThreads and CPUThreads 
apply_Htarget.h         -- Data structures -> CIJ, LIndex & RIndex Patch
CMakeLists.txt          -- CMake dependency file
matrix.cpp              -- Interface to the C code -- call Dense/Sparse
matrix.h                -- Definition for the A and B list of matrices


================================================================================
Prerequisites : 
================================================================================
cmake (version 3.7.2)
C++11 


================================================================================
Building and Running : 
================================================================================

Copy any single file from workingOpenMP / workingKokkos / workingHabaneroCPP 
to applyHTarget.cpp. 

CMakeLists.txt uses applyHTarget.cpp as the primary source file


cd miniapp_ronnie/build
cmake ../
make 

./simulateFake SYSTEM_SIZE NUM_STATES SWEEP_LOCATION NUM_K_MAT


Use the following parameters:: 

Small Test  -- System size 32  num_states = 1000  sweep location 16   k = 4

Medium Test -- System size 64   num_states = 5000   sweep location 32   k = 4

Large Test  -- System size 144   num_states = 10000  sweep location 72   k = 1



Not being used currently: 
./test_traverse CIJ_DIR X_VEC Y_VEC LIndex_Patch RIndex_Patch Sparsity_Threshold
./simulateFake CIJ_SHAPE_FILE NUM_K_MAT

================================================================================
Author : 
================================================================================
Arghya Chatterjee @chatterjeea


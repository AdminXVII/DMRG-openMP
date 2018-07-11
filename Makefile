BLAS_DIR = $(HOME)/Projects/openblas/OpenBLAS
BLAS_LIBS = -L$(BLAS_DIR) -lopenblas

KRON_LIB_DIR = util
KRON_LIBS = -L$(KRON_LIB_DIR) -lkron

OMP_FLAGS = -fopenmp

CXX = g++

TARGETS = traverse_seq simulate_driver_seq simulate_driver_omp traverse_omp \
		  traverse_indep simulate_driver_indep simulate_driver_generate
all:$(TARGETS)

working_Htarget.o : working_Htarget.cpp
	$(CXX) -c $(OMP_FLAGS) -I. -I$(KRON_LIB_DIR) $<

3Level_Indi_aHT.o : 3Level_Indi_aHT.cpp
	$(CXX) -c $(OMP_FLAGS) -I. -I$(KRON_LIB_DIR) $<

%.o : %.cpp
	$(CXX) -c -I. -I$(KRON_LIB_DIR) -I../matlab/Analysis/ $<

simulate_driver_generate:  simulate_driver_generate.o Matrix.o seq_apply_Htarget.o
	$(CXX) $(OMP_FLAGS) ../matlab/Analysis/gen_patches_comb.c $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

simulate_driver_seq:  simulate_driver.o Matrix.o seq_apply_Htarget.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

traverse_seq:  traverse.o Matrix.o seq_apply_Htarget.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

simulate_driver_omp:  simulate_driver.o Matrix.o working_Htarget.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

traverse_omp:  traverse.o Matrix.o working_Htarget.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

simulate_driver_indep:  simulate_driver.o Matrix.o 3Level_Indi_aHT.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

traverse_indep:  traverse.o Matrix.o 3Level_Indi_aHT.o
	$(CXX) $(OMP_FLAGS) $^ -o $@ $(KRON_LIBS) $(BLAS_LIBS)

clean:
	$(RM) $(TARGETS) *.o
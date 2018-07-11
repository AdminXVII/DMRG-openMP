MKL_DIR=/opt/intel/compilers_and_libraries_2017.2.174/linux/mkl/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_DIR/lib/intel64_lin
gcc -g -fopenmp   -DUSE_INTEL_MKL \
    -o test_vbatch_mkl \
    cal_kron_flops.c \
    gen_patches_comb.c \
    setup_vbatch.c \
    dmrg_malloc.c \
    dmrg_vbatch.c \
    apply_Htarget_vbatch.c  \
    test_vbatch.c \
    -I$MKL_DIR/include \
    -L$MKL_DIR/lib/intel64_lin  -lmkl_avx2 -lmkl_core -lmkl_rt -lmkl_gnu_thread  \
    -lpthread -lm


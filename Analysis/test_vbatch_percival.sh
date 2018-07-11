MKL_DIR=/opt/intel/compilers_and_libraries_2017.2.174/linux/mkl/
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:$MKL_DIR/lib/intel64_lin
icc -fast   -qopenmp -xmic-avx512 -UUSE_INTEL_MKL \
    -o test_vbatch_omp \
    cal_kron_flops.c \
    gen_patches_comb.c \
    setup_vbatch.c \
    dmrg_malloc.c \
    dmrg_vbatch.c \
    apply_Htarget_vbatch.c  \
    test_vbatch.c \
    -mkl \
    -lm 

icc -fast   -qopenmp -xmic-avx512 -DUSE_INTEL_MKL \
    -o test_vbatch_mkl \
    cal_kron_flops.c \
    gen_patches_comb.c \
    setup_vbatch.c \
    dmrg_malloc.c \
    dmrg_vbatch.c \
    apply_Htarget_vbatch.c  \
    test_vbatch.c \
    -mkl  \
    -lm 

cp test_vbatch_mkl test_vbatch_omp $MEMBERWORK/stf006/bin   

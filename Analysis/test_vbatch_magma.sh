MAGMA_DIR=/home/efdazedo/MAGMA/magma-2.2.0
# CUDA_DIR=/opt/pgi/linux86-64/2016/cuda/8.0/
CUDA_DIR=/usr/local/cuda-8.0
# CUDA_DIR=/opt/pgi/linux86-64/2017/cuda/7.5-pgprof/
g++ -g   -fopenmp -DUSE_MAGMA \
    -o test_vbatch_magma \
    cal_kron_flops.c \
    gen_patches_comb.c \
    setup_vbatch.c \
    dmrg_malloc.c \
    dmrg_vbatch.c \
    apply_Htarget_vbatch.c  \
    test_vbatch.c \
    -I$MAGMA_DIR/include \
    -L$MAGMA_DIR/lib -lmagma \
    -Wl,-rpath,$MAGMA_DIR/lib \
    -I$CUDA_DIR/include \
    -L$CUDA_DIR/lib64 -lcublas -lcusparse -lcudart \
    -Wl,-rpath,$CUDA_DIR/lib \
    -L/usr/lib -llapack \
    -L/usr/lib/openblas-base -lopenblas \
    -lm -lpthread


# pgcc -g -Mbounds -mp -Minfo=all \
touch test_vbatch_pgi
rm test_vbatch_pgi

pgcc -fast  -mp -Minfo=all \
    -o test_vbatch_pgi \
    cal_kron_flops.c \
    gen_patches_comb.c \
    setup_vbatch.c \
    dmrg_malloc.c \
    dmrg_vbatch.c \
    apply_Htarget_vbatch.c  \
    test_vbatch.c \
    -llapack -lblas -pgf77libs -pgf90libs -lm  


gcc -g -Wall -o test_gen_patches_comb -g \
    gen_patches_comb.c \
    test_gen_patches_comb.c -lm
./test_gen_patches_comb > out
octave --no-wind < test_gen_patches_comb.m > out2
diff out out2

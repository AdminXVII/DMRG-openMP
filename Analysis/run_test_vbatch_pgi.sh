source test_vbatch_pgi.sh

export OMP_STACKSIZE=2G

export OMP_NUM_THREADS=1
./test_vbatch_pgi -n 8 -m 512 -o 1   >/tmp/out.1.1

export OMP_NUM_THREADS=2
./test_vbatch_pgi -n 8 -m 512 -o 1   >/tmp/out.1.2

export OMP_NUM_THREADS=4
./test_vbatch_pgi -n 8 -m 512 -o 1   >/tmp/out.1.4



grep sd_ /tmp/out.?.?

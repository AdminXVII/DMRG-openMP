#!/bin/bash
gcc -c -O3 *.c
ar rcs lib.a *.o

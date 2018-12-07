#!/bin/bash
gcc -c -O3 *.c
ar rcs libAnalysis.a *.o

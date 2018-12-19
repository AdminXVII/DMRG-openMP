#!/bin/bash
armclang -c -O3 *.c
ar rcs lib.a *.o

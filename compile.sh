#!/bin/bash

gcc -O3 -fopenmp -fpic -c vdos.c
gcc -shared -lgomp vdos.o -o vdos.so

#!/bin/bash

gcc -O3 -fopenmp -fpic -lgsl -lgslcblas -c vdos.c
gcc -shared -lgomp -lgsl -lgslcblas vdos.o -o vdos.so

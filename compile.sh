#!/bin/bash

gcc -O3 -fpic -c vdos.c
gcc -shared vdos.o -o vdos.so

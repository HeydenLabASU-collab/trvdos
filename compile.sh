#!/bin/bash

gcc -fpic -c vdos.c
gcc -shared vdos.o -o vdos.so

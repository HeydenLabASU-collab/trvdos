#!/bin/bash

gcc -shared vdos.o \
    -lgomp \
    -lgsl \
    -lgslcblas \
    -lm \
    -o libvdos.so
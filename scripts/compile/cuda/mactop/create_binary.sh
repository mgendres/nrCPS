#!/bin/bash

SRCDIR=../../
BUILDDIR=../../
CUFFTDIR=/usr/local/cuda
CUDADIR=/usr/local/cuda

CXX=mpicxx

LDFLAGS="-L$CUFFTDIR/lib -L$CUDADIR/lib"
CXXFLAGS="-O2 -Wall -DMPICH_IGNORE_CXX_SEEK -m32"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$CUFFTDIR/include -I$CUDADIR/include"

MAIN=$1
$CXX -lcufft -lcudart $MAIN $CXXFLAGS $BUILDDIR/nrCPS.a $LDFLAGS $INCLUDE_FLAGS

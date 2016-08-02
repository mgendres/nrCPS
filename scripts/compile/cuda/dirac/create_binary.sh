#!/bin/bash

SRCDIR=../../
BUILDDIR=../../
CUFFTDIR=/usr/common/usg/cuda/3.2
CUDADIR=/usr/common/usg/cuda/3.2

CXX=mpicxx

LDFLAGS="-L$CUFFTDIR/lib64 -L$CUDADIR/lib64"
CXXFLAGS="-O2 -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$CUFFTDIR/include -I$CUDADIR/include"

MAIN=$1
$CXX -lcufft -lcudart  $MAIN $CXXFLAGS $BUILDDIR/nrCPS.a $LDFLAGS $INCLUDE_FLAGS

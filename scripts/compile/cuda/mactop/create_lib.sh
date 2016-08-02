#!/bin/bash

SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$SRCDIR/src
BUILDDIR2=$BUILDDIR/objs
CUFFTDIR=/usr/local/cuda
CUDADIR=/usr/local/cuda

CXX=mpicxx
AR=ar

LDFLAGS="-L$CUFFTDIR/lib -L$CUDADIR/lib"
CXXFLAGS="-O2 -w -DMPICH_IGNORE_CXX_SEEK -m32"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$CUFFTDIR/include -I$CUDADIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*/*.C 
nvcc -O2 $INCLUDE_FLAGS $SRCDIR2/*/cuda/*.cu -c -lcufft
$AR $ARFLAGS $BUILDDIR/nrCPS.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/nrCPS.a

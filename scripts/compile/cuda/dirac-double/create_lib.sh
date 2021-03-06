#!/bin/bash

SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$SRCDIR/src
BUILDDIR2=$BUILDDIR/objs
CUFFTDIR=/usr/common/usg/cuda/3.2
CUDADIR=/usr/common/usg/cuda/3.2

CXX=mpicxx
AR=ar

LDFLAGS="-L$CUFFTDIR/lib64 -L$CUDADIR/lib64"
CXXFLAGS="-O2 -w -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$CUFFTDIR/include -I$CUDADIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*/*.C 
nvcc -O2 $INCLUDE_FLAGS $SRCDIR2/*/cuda/*.cu -c -lcufft -arch sm_13
$AR $ARFLAGS $BUILDDIR/nrCPS.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/nrCPS.a

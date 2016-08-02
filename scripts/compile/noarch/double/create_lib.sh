#!/bin/bash

SRCDIR=$PWD
BUILDDIR=$PWD
SRCDIR2=$SRCDIR/src
BUILDDIR2=$BUILDDIR/objs
#FFTWDIR=$HOME/libraries/dfftw-3.2.2
FFTWDIR=$HOME/Documents/libraries/fftw-3.3.3


CC=mpicc
CXX=mpicxx
AR=ar

LDFLAGS="-L$FFTWDIR/lib"
CXXFLAGS="-O2 -w -DMPICH_IGNORE_CXX_SEEK"
ARFLAGS="ruv"

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

if [ -d $BUILDDIR2 ]; then
  rm -r $BUILDDIR2
fi

mkdir $BUILDDIR2
cd $BUILDDIR2

$CXX $CXXFLAGS $LDFLAGS $INCLUDE_FLAGS -c $SRCDIR2/*/*.C $SRCDIR2/*/noarch/*.C
$AR $ARFLAGS $BUILDDIR/nrCPS.a $BUILDDIR2/*.o
ranlib  $BUILDDIR/nrCPS.a

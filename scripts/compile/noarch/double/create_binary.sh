#!/bin/bash

SRCDIR=../../
BUILDDIR=../../
FFTWDIR=$HOME/libraries/dfftw-3.2.2
FFTWDIR=$HOME/Documents/libraries/fftw-3.3.3
#LAPACKDIR=$HOME/libraries/lapack

CXX=mpicxx

LDFLAGS="-L$FFTWDIR/lib -lm  -Wl,-framework -Wl,Accelerate"
CXXFLAGS="-O2 -Wall -DMPICH_IGNORE_CXX_SEEK -arch x86_64"
ARFLAGS=

INCLUDE_FLAGS="-I$SRCDIR/include -I$FFTWDIR/include"

MAIN=$1
$CXX $CXXFLAGS $MAIN $BUILDDIR/nrCPS.a $LDFLAGS $INCLUDE_FLAGS -lfftw3

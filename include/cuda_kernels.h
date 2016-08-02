#include <cuda.h>
#include <cuComplex.h>

#ifdef USE_SINGLE
typedef cuComplex Complex;
#endif
#ifdef USE_DOUBLE
typedef cuDoubleComplex Complex;
#endif
#ifdef USE_LONG_DOUBLE
#error "USE_GPU not supported with USE_LONG_DOUBLE"
#endif

#ifndef INCLUDED_CUDA_KERNELS
#define INCLUDED_CUDA_KERNELS

__global__ void Complex_Eq_cFloat(Complex*, Float*, int);
__global__ void Complex_TimesEq_Float(Complex*,  Float*, int);
__global__ void Complex_PlusEq_One_Div_Vol(Complex*, int);

__global__ void Complex_DivEq_Float(Complex*, Float* , int);
__global__ void Complex_TimesEq_Float(Complex*, Float* , int);
__global__ void Complex_TimesEq_Complex(Complex*, Complex*,  int);

__global__ void two_body_run(Float*, Complex*, Complex*, int*, Complex*, int);
__global__ void two_body_run(Float*, Complex**, Complex**, int*, Complex*, int, int);

__global__ void one_body_project(Complex*, Complex*, Complex*, int);
__global__ void inner_product(Complex**, Complex**, Complex*, int, int);

#endif

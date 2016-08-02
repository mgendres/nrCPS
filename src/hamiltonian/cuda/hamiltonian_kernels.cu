#include "config.h"
#include "cuda_kernels.h"

__global__ void Complex_Eq_cFloat(Complex* A, Float* B, int len)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < len) {
    A[tid].x = B[2*tid];
    A[tid].y = B[2*tid+1];
    tid += blockDim.x * gridDim.x;
  }

}

__global__ void Complex_DivEq_Float(Complex* A, Float* B, int len)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < len) {
    A[tid].x /= B[tid];
    A[tid].y /= B[tid];
    tid += blockDim.x * gridDim.x;
  }

}


__global__ void Complex_TimesEq_Float(Complex* A,  Float* B, int len)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < len) {
    A[tid].x *= B[tid];
    A[tid].y *= B[tid];
    tid += blockDim.x * gridDim.x;
  }

}

__global__ void Complex_TimesEq_Complex(Complex* A, Complex* B,  int len)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  Complex z;
  while (tid < len) {
    z.x = A[tid].x * B[tid].x - A[tid].y * B[tid].y;
    z.y = A[tid].y * B[tid].x + A[tid].x * B[tid].y;
    A[tid] = z;
    tid += blockDim.x * gridDim.x;
  }

}

__global__ void Complex_PlusEq_One_Div_Vol(Complex* A, int len)
{

  int tid = threadIdx.x + blockIdx.x * blockDim.x;

  while (tid < len) {
    A[tid].x += 1.0;
    A[tid].x /= len;
    A[tid].y /= len;
    tid += blockDim.x * gridDim.x;
  }

}



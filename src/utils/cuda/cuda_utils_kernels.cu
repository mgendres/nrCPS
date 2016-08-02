#include "config.h"
#include "cuda_kernels.h"

__global__ void  inner_product(Complex** kets, Complex** bras, Complex* matrix, int nf, int vol)
{

  __shared__ Complex cache[THREADS];

  int col = blockIdx.x; 
  int row = blockIdx.y; 

  int tid=threadIdx.x;
  int cid=threadIdx.x;

  int i;

  cache[cid].x = 0.0;
  cache[cid].y = 0.0;

  Complex* ket = kets[row];
  Complex* bra = bras[col];

  while( tid < vol) {

    cache[cid].x += ket[tid].x * bra[tid].x + ket[tid].y * bra[tid].y ;
    cache[cid].y += ket[tid].x * bra[tid].y - ket[tid].y * bra[tid].x ;

    tid += blockDim.x;
  }

  __syncthreads();

  i = blockDim.x/2;

  while (i != 0) {
    if (cid < i) {
      cache[cid].x += cache[cid + i].x;
      cache[cid].y += cache[cid + i].y;
    }
    __syncthreads();
    i /= 2;
  }

  if (cid == 0) {
    matrix[col+row*nf] = cache[0];
  }

}


#include "cuda_utils.h"
#include "cuda_kernels.h"

__global__ void one_body_project(Complex* psi, Complex* phi, Complex* psum, int len)
{

  __shared__ Complex cache[THREADS];

  int tid = threadIdx.x + blockIdx.x * blockDim.x; // thread id
  int cid = threadIdx.x; // cache id

  cache[cid].x = 0.0;
  cache[cid].y = 0.0;


  while (tid < len) {
    cache[cid].x += (psi[tid].x * phi[tid].x + psi[tid].y * phi[tid].y );
    cache[cid].y += (psi[tid].x * phi[tid].y - psi[tid].y * phi[tid].x );
    tid += blockDim.x * gridDim.x;    
  }

  __syncthreads();

  int i = blockDim.x/2;

  while (i != 0) {
    if (cid < i) {
      cache[cid].x += cache[cid + i].x;
      cache[cid].y += cache[cid + i].y;
    }
    __syncthreads();
    i /= 2;
  }

  if (cid == 0) {
    psum[blockIdx.x] = cache[0];
  }

}

#include "config.h"
#include "cuda_kernels.h"

__global__ void two_body_run(Float* wavefunc, Complex* prop1, Complex* prop2, int* opposite_parity_index, Complex* psum, int vol)
{

  __shared__ Complex cache[THREADS];

  int tid = threadIdx.x + blockIdx.x * blockDim.x; // thread id
  int cid = threadIdx.x; // cache id

  cache[cid].x = 0.0;
  cache[cid].y = 0.0;

  int k;

  while (tid < vol) {
    k = opposite_parity_index[tid];
    cache[cid].x += (prop1[tid].x * prop2[k].x - prop1[tid].y * prop2[k].y ) * wavefunc[tid];
    cache[cid].y += (prop1[tid].x * prop2[k].y + prop1[tid].y * prop2[k].x ) * wavefunc[tid];
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


/*
__global__ void  two_body_run(Float* wavefunc, Complex** props1, Complex** props2, int* opposite_parity_index, Complex* slater, int nf, int vol)
{

  __shared__ Complex cache[THREADS];

  int col = blockIdx.x; 
  int row = blockIdx.y; 

  int tid=threadIdx.x;
  int cid=threadIdx.x;

  int i;

  cache[cid].x = 0.0;
  cache[cid].y = 0.0;

  Complex* p1 = props1[row];
  Complex* p2 = props2[col];

  while( tid < vol) {

    i = opposite_parity_index[tid];
    cache[cid].x += (p1[tid].x * p2[i].x - p1[tid].y * p2[i].y ) * wavefunc[tid];
    cache[cid].y += (p1[tid].x * p2[i].y + p1[tid].y * p2[i].x ) * wavefunc[tid];

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
    slater[col+row*nf] = cache[0];
  }

}
*/



/*
#define Q 2
__global__ void  two_body_run(Float* wavefunc, Complex** props1, Complex** props2, int* opposite_parity_index, Complex* slater, int nf, int vol)
{

  __shared__ Complex cache[Q][Q][THREADS];

  int col = Q*blockIdx.x; 
  int row = Q*blockIdx.y; 

  int tid=threadIdx.x;
  int cid=threadIdx.x;

  int i;

  for (int s=0; s<Q; ++s) {
    for (int t=0; t<Q; ++t) {
      cache[s][t][cid].x = 0.0;
      cache[s][t][cid].y = 0.0;
    }
  }


  Complex p1[Q];
  Complex p2[Q];
  Float wf;

  while( tid < vol) {

    i = opposite_parity_index[tid];

    wf = wavefunc[tid];

    for (int s=0; s<Q; ++s) {
      p1[s] = props1[row+s][tid];
    }

    for (int t=0; t<Q; ++t) {
      p2[t] =  props2[col+t][i];
    }

    for (int s=0; s<Q; ++s) {
      for (int t=0; t<Q; ++t) {
        cache[s][t][cid].x += (p1[s].x * p2[t].x - p1[s].y * p2[t].y ) * wf;
        cache[s][t][cid].y += (p1[s].x * p2[t].y + p1[s].y * p2[t].x ) * wf;
      }
    }
    tid += blockDim.x;
  }

  __syncthreads();

  i = blockDim.x/2;

  while (i != 0) {
    if (cid < i) {

      for (int s=0; s<Q; ++s) {
        for (int t=0; t<Q; ++t) {
          cache[s][t][cid].x += cache[s][t][cid + i].x;
          cache[s][t][cid].y += cache[s][t][cid + i].y;
        }
      }

    }
    __syncthreads();
    i /= 2;
  }

  if (cid == 0) {

    for (int s=0; s<Q; ++s) {
      for (int t=0; t<Q; ++t) {
        slater[(col+t)+(row+s)*nf] = cache[s][t][0];
      }
    }

  }

}
*/




/*
#define QR 3
#define QC 2
__global__ void  two_body_run(Float* wavefunc, Complex** props1, Complex** props2, int* opposite_parity_index, Complex* slater, int nf, int vol)
{

  __shared__ Complex cache[THREADS][QR][QC];

  int col = QR*blockIdx.x; 
  int row = QC*blockIdx.y; 

  int tid=threadIdx.x;
  int cid=threadIdx.x;

  int i;

  for (int s=0; s<QR; ++s) {
    for (int t=0; t<QC; ++t) {
      cache[cid][s][t].x = 0.0;
      cache[cid][s][t].y = 0.0;
    }
  }


  Complex p1[QR];
  Complex p2[QC];
  Float wf;

  while( tid < vol) {

    i = opposite_parity_index[tid];

    wf = wavefunc[tid];

    for (int s=0; s<QR; ++s) {
      p1[s] = props1[row+s][tid];
    }

    for (int t=0; t<QC; ++t) {
      p2[t] =  props2[col+t][i];
    }

    for (int s=0; s<QR; ++s) {
      //for (int t=0; t<QC; ++t) {
        cache[cid][s][0].x += (p1[s].x * p2[0].x - p1[s].y * p2[0].y ) * wf;
        cache[cid][s][0].y += (p1[s].x * p2[0].y + p1[s].y * p2[0].x ) * wf;
        cache[cid][s][1].x += (p1[s].x * p2[1].x - p1[s].y * p2[1].y ) * wf;
        cache[cid][s][1].y += (p1[s].x * p2[1].y + p1[s].y * p2[1].x ) * wf;
      //}
    }

    tid += blockDim.x;

  }

  __syncthreads();

  i = blockDim.x/2;

  while (i != 0) {
    if (cid < i) {

      for (int s=0; s<QR; ++s) {
        //for (int t=0; t<QC; ++t) {
          cache[cid][s][0].x += cache[cid+i][s][0].x;
          cache[cid][s][0].y += cache[cid+i][s][0].y;
          cache[cid][s][1].x += cache[cid+i][s][1].x;
          cache[cid][s][1].y += cache[cid+i][s][1].y;
        //}
      }

    }
    __syncthreads();
    i /= 2;
  }

  if (cid == 0) {

    for (int s=0; s<QR; ++s) {
      //for (int t=0; t<QC; ++t) {
        slater[(col+0)+(row+s)*nf] = cache[0][s][0];
        slater[(col+1)+(row+s)*nf] = cache[0][s][1];
      //}
    }

  }

}
*/




#define QR 3
#define QC 3
__global__ void  two_body_run(Float* wavefunc, Complex** props1, Complex** props2, int* opposite_parity_index, Complex* slater, int nf, int vol)
{

  __shared__ Complex cache[QR][QC][THREADS];

  int col = QR*blockIdx.x; 
  int row = QC*blockIdx.y; 

  int tid=threadIdx.x;
  int cid=threadIdx.x;

  int i;

  cache[0][0][cid].x = 0.0;
  cache[0][0][cid].y = 0.0;
  cache[0][1][cid].x = 0.0;
  cache[0][1][cid].y = 0.0;
  cache[0][2][cid].x = 0.0;
  cache[0][2][cid].y = 0.0;
  cache[1][0][cid].x = 0.0;
  cache[1][0][cid].y = 0.0;
  cache[1][1][cid].x = 0.0;
  cache[1][1][cid].y = 0.0;
  cache[1][2][cid].x = 0.0;
  cache[1][2][cid].y = 0.0;
  cache[2][0][cid].x = 0.0;
  cache[2][0][cid].y = 0.0;
  cache[2][1][cid].x = 0.0;
  cache[2][1][cid].y = 0.0;
  cache[2][2][cid].x = 0.0;
  cache[2][2][cid].y = 0.0;

  Complex p1[QR];
  Complex p2[QC];
  Float wf;

  while( tid < vol) {

    i = opposite_parity_index[tid];

    wf = wavefunc[tid];

    p1[0] = props1[row+0][tid];
    p1[1] = props1[row+1][tid];
    p1[2] = props1[row+2][tid];

    p2[0] =  props2[col+0][i];
    p2[1] =  props2[col+1][i];
    p2[2] =  props2[col+2][i];

//  if (row+0<nf) { p1[0] = props1[row+0][tid]; } else { p1[0].x=0.0; p1[0].y=0.0; }
//  if (row+1<nf) { p1[1] = props1[row+1][tid]; } else { p1[1].x=0.0; p1[1].y=0.0; }
//  if (row+2<nf) { p1[2] = props1[row+2][tid]; } else { p1[2].x=0.0; p1[2].y=0.0; }

//  if (col+0<nf) { p2[0] = props2[col+0][i]; } else { p2[0].x=0.0; p2[0].y=0.0; }
//  if (col+1<nf) { p2[1] = props2[col+1][i]; } else { p2[1].x=0.0; p2[1].y=0.0; }
//  if (col+2<nf) { p2[2] = props2[col+2][i]; } else { p2[2].x=0.0; p2[2].y=0.0; }

    cache[0][0][cid].x += (p1[0].x * p2[0].x - p1[0].y * p2[0].y ) * wf;
    cache[0][0][cid].y += (p1[0].x * p2[0].y + p1[0].y * p2[0].x ) * wf;
    cache[0][1][cid].x += (p1[0].x * p2[1].x - p1[0].y * p2[1].y ) * wf;
    cache[0][1][cid].y += (p1[0].x * p2[1].y + p1[0].y * p2[1].x ) * wf;
    cache[0][2][cid].x += (p1[0].x * p2[2].x - p1[0].y * p2[2].y ) * wf;
    cache[0][2][cid].y += (p1[0].x * p2[2].y + p1[0].y * p2[2].x ) * wf;
    cache[1][0][cid].x += (p1[1].x * p2[0].x - p1[1].y * p2[0].y ) * wf;
    cache[1][0][cid].y += (p1[1].x * p2[0].y + p1[1].y * p2[0].x ) * wf;
    cache[1][1][cid].x += (p1[1].x * p2[1].x - p1[1].y * p2[1].y ) * wf;
    cache[1][1][cid].y += (p1[1].x * p2[1].y + p1[1].y * p2[1].x ) * wf;
    cache[1][2][cid].x += (p1[1].x * p2[2].x - p1[1].y * p2[2].y ) * wf;
    cache[1][2][cid].y += (p1[1].x * p2[2].y + p1[1].y * p2[2].x ) * wf;
    cache[2][0][cid].x += (p1[2].x * p2[0].x - p1[2].y * p2[0].y ) * wf;
    cache[2][0][cid].y += (p1[2].x * p2[0].y + p1[2].y * p2[0].x ) * wf;
    cache[2][1][cid].x += (p1[2].x * p2[1].x - p1[2].y * p2[1].y ) * wf;
    cache[2][1][cid].y += (p1[2].x * p2[1].y + p1[2].y * p2[1].x ) * wf;
    cache[2][2][cid].x += (p1[2].x * p2[2].x - p1[2].y * p2[2].y ) * wf;
    cache[2][2][cid].y += (p1[2].x * p2[2].y + p1[2].y * p2[2].x ) * wf;

    tid += blockDim.x;

  }

  __syncthreads();

  i = blockDim.x/2;

  while (i != 0) {
    if (cid < i) {

      cache[0][0][cid].x += cache[0][0][cid+i].x;
      cache[0][0][cid].y += cache[0][0][cid+i].y;
      cache[0][1][cid].x += cache[0][1][cid+i].x;
      cache[0][1][cid].y += cache[0][1][cid+i].y;
      cache[0][2][cid].x += cache[0][2][cid+i].x;
      cache[0][2][cid].y += cache[0][2][cid+i].y;
      cache[1][0][cid].x += cache[1][0][cid+i].x;
      cache[1][0][cid].y += cache[1][0][cid+i].y;
      cache[1][1][cid].x += cache[1][1][cid+i].x;
      cache[1][1][cid].y += cache[1][1][cid+i].y;
      cache[1][2][cid].x += cache[1][2][cid+i].x;
      cache[1][2][cid].y += cache[1][2][cid+i].y;
      cache[2][0][cid].x += cache[2][0][cid+i].x;
      cache[2][0][cid].y += cache[2][0][cid+i].y;
      cache[2][1][cid].x += cache[2][1][cid+i].x;
      cache[2][1][cid].y += cache[2][1][cid+i].y;
      cache[2][2][cid].x += cache[2][2][cid+i].x;
      cache[2][2][cid].y += cache[2][2][cid+i].y;

    }
    __syncthreads();
    i /= 2;
  }

  if ( (THREADS==96)||(THREADS==112) ) {

    cache[0][0][0].x += cache[0][0][2].x;
    cache[0][0][0].y += cache[0][0][2].y;
    cache[0][1][0].x += cache[0][1][2].x;
    cache[0][1][0].y += cache[0][1][2].y;
    cache[0][2][0].x += cache[0][2][2].x;
    cache[0][2][0].y += cache[0][2][2].y;
    cache[1][0][0].x += cache[1][0][2].x;
    cache[1][0][0].y += cache[1][0][2].y;
    cache[1][1][0].x += cache[1][1][2].x;
    cache[1][1][0].y += cache[1][1][2].y;
    cache[1][2][0].x += cache[1][2][2].x;
    cache[1][2][0].y += cache[1][2][2].y;
    cache[2][0][0].x += cache[2][0][2].x;
    cache[2][0][0].y += cache[2][0][2].y;
    cache[2][1][0].x += cache[2][1][2].x;
    cache[2][1][0].y += cache[2][1][2].y;
    cache[2][2][0].x += cache[2][2][2].x;
    cache[2][2][0].y += cache[2][2][2].y;

  }
  __syncthreads();

  if ( THREADS==112 ) {

    cache[0][0][0].x += cache[0][0][6].x;
    cache[0][0][0].y += cache[0][0][6].y;
    cache[0][1][0].x += cache[0][1][6].x;
    cache[0][1][0].y += cache[0][1][6].y;
    cache[0][2][0].x += cache[0][2][6].x;
    cache[0][2][0].y += cache[0][2][6].y;
    cache[1][0][0].x += cache[1][0][6].x;
    cache[1][0][0].y += cache[1][0][6].y;
    cache[1][1][0].x += cache[1][1][6].x;
    cache[1][1][0].y += cache[1][1][6].y;
    cache[1][2][0].x += cache[1][2][6].x;
    cache[1][2][0].y += cache[1][2][6].y;
    cache[2][0][0].x += cache[2][0][6].x;
    cache[2][0][0].y += cache[2][0][6].y;
    cache[2][1][0].x += cache[2][1][6].x;
    cache[2][1][0].y += cache[2][1][6].y;
    cache[2][2][0].x += cache[2][2][6].x;
    cache[2][2][0].y += cache[2][2][6].y;

  }
  __syncthreads();

  if (cid == 0) {

    slater[(col+0)+(row+0)*nf] = cache[0][0][0];
    slater[(col+1)+(row+0)*nf] = cache[0][1][0];
    slater[(col+2)+(row+0)*nf] = cache[0][2][0];
    slater[(col+0)+(row+1)*nf] = cache[1][0][0];
    slater[(col+1)+(row+1)*nf] = cache[1][1][0];
    slater[(col+2)+(row+1)*nf] = cache[1][2][0];
    slater[(col+0)+(row+2)*nf] = cache[2][0][0];
    slater[(col+1)+(row+2)*nf] = cache[2][1][0];
    slater[(col+2)+(row+2)*nf] = cache[2][2][0];

//  if (row+0<nf) {
//    if (col+0<nf) slater[(col+0)+(row+0)*nf] = cache[0][0][0];
//    if (col+1<nf) slater[(col+1)+(row+0)*nf] = cache[0][0][1];
//    if (col+2<nf) slater[(col+2)+(row+0)*nf] = cache[0][0][2];
//  }

//  if (row+1<nf) {
//    if (col+0<nf) slater[(col+0)+(row+1)*nf] = cache[0][1][0];
//    if (col+1<nf) slater[(col+1)+(row+1)*nf] = cache[0][1][1];
//    if (col+2<nf) slater[(col+2)+(row+1)*nf] = cache[0][1][2];
//  }

//  if (row+2<nf) {
//    if (col+0<nf) slater[(col+0)+(row+2)*nf] = cache[0][2][0];
//    if (col+1<nf) slater[(col+1)+(row+2)*nf] = cache[0][2][1];
//    if (col+2<nf) slater[(col+2)+(row+2)*nf] = cache[0][2][2];
//  }

  }

}




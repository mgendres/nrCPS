#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "arg.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"
#include "dispersion.h"
#include "special_functions.h"

#include "cuda_utils.h"
#include "cuda_kernels.h"

complex<Float> TwoBody::Run(Float* prop1, Float* prop2)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*)";
//  VRB.Func(fname);

  two_body_run<<<BLOCKS,THREADS>>>(dev_wavefunc, (Complex*) prop1, (Complex*) prop2, dev_opposite_parity_index, (Complex*) dev_psum, vol);

  Cuda::MemCopy(psum, dev_psum, 2*BLOCKS*sizeof(Float), cudaMemcpyDeviceToHost);

  Float re=0.0;
  Float im=0.0;

  for (int i=0; i<2*BLOCKS; i+=2) {
    re += psum[i];
    im += psum[i+1];
  }

  return complex<Float>(re,im);

}


void TwoBody::Run(Float** props1, Float** props2, complex<Float>* slater , int nf)
{

  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*, complex<Float>, int)";
//  VRB.Func(fname);

  Complex* host_slater; // Eventualy move this to constructor?
  host_slater = (Complex *) malloc(nf*nf*sizeof(Complex) );

  Complex* dev_slater; // Eventualy move this to constructor?
  Cuda::Malloc( (void**)&dev_slater, nf*nf*sizeof(Complex) );

  int QR=3;
  int QC=3;
  //dim3 blocks( (nf+QR-1)/QR, (nf+QC-1)/QC);
  dim3 blocks( nf/QR, nf/QC );
  dim3 threads(THREADS,1);

  float elapsed_time;
  cudaEvent_t start, stop;
  cudaEventCreate(&start);
  cudaEventCreate(&stop);
  cudaEventRecord(start, 0);

  two_body_run<<<blocks,threads>>>(dev_wavefunc, (Complex**) props1, (Complex**) props2, dev_opposite_parity_index, dev_slater, nf, vol);
  
  cudaEventRecord(stop, 0);
  cudaEventSynchronize(stop);
  cudaEventElapsedTime(&elapsed_time, start, stop) ;
  cudaEventDestroy(start);
  cudaEventDestroy(stop);
  VRB.Flow(fname, "Total elapsed time is: %3.1f ms.", elapsed_time);
  VRB.Flow(fname, "Total bandwidth is: %3.1f GB/s.", ( (nf/QR)*(nf/QC)*GJP.Vol()*(2*QC+2*QR+1)*sizeof(double) )/( 1000000.0*elapsed_time ));


  Cuda::MemCopy(host_slater, dev_slater, nf*nf*sizeof(Complex), cudaMemcpyDeviceToHost);

  for (int i=0; i<nf*nf; ++i) { slater[i] = complex<Float>( host_slater[i].x, host_slater[i].y); }

  Cuda::Free(dev_slater);
  free(host_slater);

}



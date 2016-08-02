#include <complex>
#include <config.h>
#include <cuda.h>
#include <driver_types.h>
#ifndef INCLUDED_CUDA_UTILS
#define INCLUDED_CUDA_UTILS

namespace Cuda {

  void Malloc(void**, size_t);
  void HostAlloc(void**, size_t , unsigned int );
  void Free(void*);
  void FreeHost(void*);
  void MemCopy(void*, void*, size_t, cudaMemcpyKind);
  void DeviceQ();


  void InnerProduct(Float**, Float**, std::complex<Float>*, int, int);

}


#endif

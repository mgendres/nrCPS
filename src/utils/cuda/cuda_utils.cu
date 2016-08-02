#include "cuda_utils.h"
#include "cuda_kernels.h"
#include <stdio.h>
#include "verbose.h"
#include "error.h"
#include <sysfunc.h>
#include <sysio.h>
using namespace std;

void Cuda::Malloc(void **devPtr, size_t size)
{
  const char* fname = "void Cuda::Malloc(void**, size_t)";
  cudaError_t err = cudaMalloc(devPtr, size );
  if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }
}

void Cuda::HostAlloc(void **pHost, size_t size, unsigned int flags )
{
  const char* fname = "void Cuda::HostMalloc(void**, size_t , unsigned int )";
  cudaError_t err = cudaHostAlloc(pHost, size, flags);
  if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }
}

void Cuda::Free(void *devPtr)
{
  const char* fname = "void Cuda::Free(void*)";
  cudaError_t err = cudaFree(devPtr);
  if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }
}

void Cuda::FreeHost(void *pHost)
{
  const char* fname = "void Cuda::FreeHost(void*)";
  cudaError_t err = cudaFreeHost(pHost);
  if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }
}

void Cuda::MemCopy(void* dst, void* src, size_t size, cudaMemcpyKind kind)
{
  const char* fname = "void Cuda::MemCopy(void*, void*, size_t, cudaMemcpyKind)";
  cudaError_t err = cudaMemcpy(dst, src, size, kind);
  if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }
}

void Cuda::DeviceQ()
{
  const char* fname = "void Cuda::DeviceQ()";
  cudaError_t err;
  cudaDeviceProp device_prop;
  int count;

  for (int i=0; i<Comms::Size(); ++i) {

    if ( i == Comms::Rank() ) {

      err = cudaGetDeviceCount(&count);
      if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }

      VRB.Flow(fname, "Number of GPU devices on host %d: %d", i, count);
 
      for (int j=0; j<count; ++j) {

        err = cudaGetDeviceProperties( &device_prop, j);
        if (cudaSuccess != err) { ERR.General(fname, "Cuda error: %s\n", cudaGetErrorString(err)); }

        VRB.Flow(fname, "Host %d/device %d: Name is %s", i, j, device_prop.name );
        VRB.Flow(fname, "Host %d/device %d: Compute capability is %d.%d", i, j, device_prop.major, device_prop.minor );
        VRB.Flow(fname, "Host %d/device %d: Clock Rate is %d", i, j, device_prop.clockRate );
        VRB.Flow(fname, "Host %d/device %d: Total global memory is %ld", i, j, device_prop.totalGlobalMem );
        VRB.Flow(fname, "Host %d/device %d: Total constant memory is %ld", i, j, device_prop.totalConstMem );
        VRB.Flow(fname, "Host %d/device %d: Multiprocessor count is %d", i, j, device_prop.multiProcessorCount );
        VRB.Flow(fname, "Host %d/device %d: Shared memory per block %ld", i, j, device_prop.sharedMemPerBlock );
        VRB.Flow(fname, "Host %d/device %d: Max threads per block is %d", i, j, device_prop.maxThreadsPerBlock );
        VRB.Flow(fname, "Host %d/device %d: Max thread dimensions is (%d, %d, %d)", i, j, device_prop.maxThreadsDim[0],
                                                                                           device_prop.maxThreadsDim[1],
                                                                                           device_prop.maxThreadsDim[2] );
        VRB.Flow(fname, "Host %d/device %d: Max grid dimensions is (%d, %d, %d)\n", i, j, device_prop.maxGridSize[0],
                                                                                         device_prop.maxGridSize[1],
                                                                                         device_prop.maxGridSize[2] );
      }

    }
    Comms::Sync();
  }

}

void Cuda::InnerProduct(Float** kets, Float** bras, complex<Float>* mat, int nf, int vol)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*, complex<Float>, int)";
//  VRB.Func(fname);

  Complex* host_mat; // Eventualy move this to constructor?
  host_mat= (Complex *) malloc(nf*nf*sizeof(Complex) );

  Complex* dev_mat; // Eventualy move this to constructor?
  Cuda::Malloc( (void**)&dev_mat, nf*nf*sizeof(Complex) );

//float elapsed_time;
//cudaEvent_t start, stop;
//cudaEventCreate(&start);
//cudaEventCreate(&stop);
//cudaEventRecord(start, 0);
 

  dim3 blocks(nf,nf);
  dim3 threads(THREADS,1);
  inner_product<<<blocks,threads>>>( (Complex**) kets, (Complex**) bras, dev_mat, nf, vol);

  
//cudaEventRecord(stop, 0);
//cudaEventSynchronize(stop);
//cudaEventElapsedTime(&elapsed_time, start, stop) ;
//printf("Total elapsed time is: %3.1f ms.\n", elapsed_time);
//cudaEventDestroy(start);
//cudaEventDestroy(stop);
  

  Cuda::MemCopy(host_mat, dev_mat, nf*nf*sizeof(Complex), cudaMemcpyDeviceToHost);

  for (int i=0; i<nf*nf; ++i) { mat[i] = complex<Float>( host_mat[i].x, host_mat[i].y); }

  Cuda::Free(dev_mat);
  free(host_mat);

}


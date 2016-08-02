#include <fourier.h>
#include "global_job_parameter.h"
#include "verbose.h"
#include "error.h"

#include "cuda_utils.h"

//---- Some file scoped variables
Fourier::Fourier(int n)
{

  const char* fname = "void Fourier::Initialize()";

  VRB.Debug(fname, "Allocating memory and creating plans for CUFFT.");

  batch_size = n;
  Cuda::Malloc((void**)&b, batch_size*GJP.Vol()*sizeof(fftComplex));

  int dims[3] = { GJP.Xsites() , GJP.Ysites(), GJP.Zsites() };
#ifdef USE_SINGLE
  cufftPlanMany(&plan, 3, dims , NULL, 1, 0, NULL, 1, 0, CUFFT_C2C, batch_size);
#endif
#ifdef USE_DOUBLE
  cufftPlanMany(&plan, 3, dims , NULL, 1, 0, NULL, 1, 0, CUFFT_Z2Z, batch_size);
#endif
#ifdef USE_LONG_DOUBLE
#error "USE_GPU not supported with USE_LONG_DOUBLE."
#endif

}

Fourier::~Fourier()
{

  const char* fname = "void Fourier::Finalize()";

  VRB.Debug(fname, "Deallocating memory and destoying plans for CUFFT.");
  Cuda::Free(b);
  cufftDestroy(plan);

}

void Fourier::Forward()
{

  const char* fname = "void Fourier::Forward()";

#ifdef USE_SINGLE
  cufft_result = cufftExecC2C(plan, b, b, CUFFT_FORWARD);
#endif
#ifdef USE_DOUBLE
  cufft_result = cufftExecZ2Z(plan, b, b, CUFFT_FORWARD);
#endif
#ifdef USE_LONG_DOUBLE
#error "USE_GPU not supported with USE_LONG_DOUBLE."
#endif

  if (cufft_result != CUFFT_SUCCESS) { ERR.General(fname, "Failed to CUFFT; error code: %d.", cufft_result); }
    
}

void Fourier::Backward()
{

  const char* fname = "void Fourier::Backward()";

#ifdef USE_SINGLE
  cufft_result = cufftExecC2C(plan, b, b, CUFFT_INVERSE);
#endif
#ifdef USE_DOUBLE
  cufft_result = cufftExecZ2Z(plan, b, b, CUFFT_INVERSE);
#endif
#ifdef USE_LONG_DOUBLE
#error "USE_GPU not supported with USE_LONG_DOUBLE."
#endif

  if (cufft_result != CUFFT_SUCCESS) { ERR.General(fname, "Failed to CUFFT; error code: %d.", cufft_result); }

}

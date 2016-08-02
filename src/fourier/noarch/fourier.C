#include <fourier.h>
#include "global_job_parameter.h"
#include "verbose.h"
#include "error.h"

Fourier::Fourier(int n)
{

  const char* fname = "void Fourier::Initialize()";

  VRB.Debug(fname, "Allocating memory and creating plans for FFTW.");

  batch_size = n;

#ifdef USE_SINGLE
  b = (fftComplex*) fftwf_malloc(batch_size*GJP.Vol()*sizeof(fftComplex));
#endif
#ifdef USE_DOUBLE
  b = (fftComplex*) fftw_malloc(batch_size*GJP.Vol()*sizeof(fftComplex));
#endif
#ifdef USE_LONG_DOUBLE
  b = (fftComplex*) fftwl_malloc(batch_size*GJP.Vol()*sizeof(fftComplex));
#endif

  // Below needs to be adjusted for batch ffts; double check in place
  int vol = GJP.Vol();
  int dims[3] = { GJP.Xsites(), GJP.Ysites(), GJP.Zsites() };

#ifdef USE_SINGLE
  p1 = fftwf_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_BACKWARD, FFTW_EXHAUSTIVE); 
  p2 = fftwf_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_FORWARD, FFTW_EXHAUSTIVE); 
#endif
#ifdef USE_DOUBLE
  p1 = fftw_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_BACKWARD, FFTW_EXHAUSTIVE); 
  p2 = fftw_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_FORWARD, FFTW_EXHAUSTIVE); 
#endif
#ifdef USE_LONG_DOUBLE
  p1 = fftwl_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_BACKWARD, FFTW_EXHAUSTIVE); 
  p2 = fftwl_plan_many_dft(3, dims, batch_size, b, NULL, 1, vol, b, NULL, 1, vol, FFTW_FORWARD, FFTW_EXHAUSTIVE); 
#endif


}

Fourier::~Fourier()
{

  const char* fname = "void Fourier::Finalize()";

    //---- This is the last instance of Hamiltonian
    VRB.Debug(fname, "Deallocating memory and destoying plans for FFTW.");

#ifdef USE_SINGLE
    fftwf_free(b);
    fftwf_destroy_plan(p1);
    fftwf_destroy_plan(p2);
    fftwf_cleanup();
#endif
#ifdef USE_DOUBLE
    fftw_free(b);
    fftw_destroy_plan(p1);
    fftw_destroy_plan(p2);
    fftw_cleanup();
#endif
#ifdef USE_LONG_DOUBLE
    fftwl_free(b);
    fftwl_destroy_plan(p1);
    fftwl_destroy_plan(p2);
    fftwl_cleanup();
#endif


}

void Fourier::Forward()
{

#ifdef USE_SINGLE
  fftwf_execute(p2);
#endif
#ifdef USE_DOUBLE
  fftw_execute(p2);
#endif
#ifdef USE_LONG_DOUBLE
  fftwl_execute(p2);
#endif

}

void Fourier::Backward()
{
#ifdef USE_SINGLE
  fftwf_execute(p1);
#endif
#ifdef USE_DOUBLE
  fftw_execute(p1);
#endif
#ifdef USE_LONG_DOUBLE
  fftwl_execute(p1);
#endif
}

#include "config.h"
#include "fourier_typedefs.h"

#ifndef INCLUDED_FOURIER
#define INCLUDED_FOURIER

class Fourier
{
  private:
    int batch_size;

    fftComplex* b;
#ifdef USE_GPU
    cufftHandle plan;
    cufftResult cufft_result;
#else

#ifdef USE_SINGLE
    fftwf_plan p1;
    fftwf_plan p2;
#endif
#ifdef USE_DOUBLE
    fftw_plan p1;
    fftw_plan p2;
#endif
#ifdef USE_LONG_DOUBLE
    fftwl_plan p1;
    fftwl_plan p2;
#endif

#endif

    Fourier& operator=(const Fourier&);
    Fourier(const Fourier&);

  public:
    explicit Fourier(int);
    ~Fourier();

    int BatchSize();
    void Forward();
    void Backward();
    fftComplex* Get(int);
    //Phase phase;
};


#endif

#ifdef USE_GPU
#include "cufft.h"
#else
#include "fftw3.h"
#endif

#if defined(USE_GPU) && defined(USE_SINGLE)
typedef cufftComplex fftComplex;
#endif

#if defined(USE_GPU) && defined(USE_DOUBLE)
typedef cufftDoubleComplex fftComplex;
#endif

#if defined(USE_GPU) && defined(USE_LONG_DOUBLE)
#error "USE_GPU not supported with USE_LONG_DOUBLE."
#endif

#if !defined(USE_GPU) && defined(USE_SINGLE)
typedef fftwf_complex fftComplex;
#endif

#if !defined(USE_GPU) && defined(USE_DOUBLE)
typedef fftw_complex fftComplex;
#endif

#if !defined(USE_GPU) && defined(USE_LONG_DOUBLE)
typedef fftwl_complex fftComplex;
#endif




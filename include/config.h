//---- Global options for compiling nrCPS

#ifndef INCLUDED_CONFIG
#define INCLUDED_CONFIG

//---- Version information (printed at beginning of simulation)
#define VERSION_MAJOR 2
#define VERSION_MINOR 7
#define VERSION_SUB 6
#define VERSION_STR "nrCPS_v2_7_6"

//---- Comment out to disable MPI support
#define USE_MPI

//---- Comment out to disable GPU support
//#define USE_GPU

//---- Uncomment desired precision level
//---- NOTE: if GPU is enabled then the precision of the host must be the same as that on the device
//---- For double precision, must compile *.cu code using -arch sm_13 option
//#define USE_SINGLE
#define USE_DOUBLE
//#define USE_LONG_DOUBLE

//-------------------------------------------------------------------------
//---- DO NOT EDIT BELOW THIS LINE UNLESS YOU KNOW WHAT YOU ARE DOING -----
//-------------------------------------------------------------------------

//----  GPU specific properties: specify the number of block BLOCKS and the number of threads THREADS
#ifdef USE_GPU
#define BLOCKS  128
#define THREADS 112 // Should be a multiple of 32 for best performance
#endif

//---- Set the precision; precision is also specified in fourier.h and cuda_kernels.h

#if defined USE_SINGLE
const int PREC=6;  // Specifies precision of output
typedef float Float;
#elif defined USE_DOUBLE
const int PREC=15;  // Specifies precision of output
typedef double Float;
#elif defined USE_LONG_DOUBLE
const int PREC=18;  // Specifies precision of output
typedef long double Float;
#else
#error "Precision not supported; specify either USE_SINGLE, USE_DOUBLE or USE_LONG_DOUBLE."
#endif

#endif

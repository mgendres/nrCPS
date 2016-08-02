#include <fourier.h>
#include "global_job_parameter.h"
#include "verbose.h"
#include "error.h"

int Fourier::BatchSize()
{
  return batch_size;
}

fftComplex* Fourier::Get(int n)
{
  const char* fname = "fftComplex* Fourier::Get(int)";
  if (n<0) { ERR.General(fname,"Argument must be non-negative."); }
  if (n>=batch_size) { ERR.General(fname,"Argument must be less than batch_size."); }
  return b + n*GJP.Vol();
}

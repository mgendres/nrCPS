#include <complex>
#include "config.h"
#include "enum.h"

#ifndef INCLUDED_SYSFUNC
#define INCLUDED_SYSFUNC

namespace Comms {

  extern bool initialized;

  void Initialize();
  void Finalize();
  void RaiseError();

  int Rank();

  int Size();

  void Sync();

  std::complex<Float> GlobalSum(std::complex<Float>);
  void GlobalSum(std::complex<Float>*, int);

}

#endif

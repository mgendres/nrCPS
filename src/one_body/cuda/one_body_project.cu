#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "arg.h"
#include "one_body.h"
#include "verbose.h"
#include "special_functions.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"

#include "cuda_utils.h"
#include "cuda_kernels.h"

complex<Float> OneBody::Project(Float* phi)
{

  one_body_project<<<BLOCKS,THREADS>>>( (Complex*) dev_psi, (Complex*) phi, (Complex*) dev_psum, vol);

  Cuda::MemCopy(psum, dev_psum, 2*BLOCKS*sizeof(Float), cudaMemcpyDeviceToHost);

  Float re=0.0;
  Float im=0.0;

  for (int i=0; i<2*BLOCKS; i+=2) {
    re += psum[i];
    im += psum[i+1];
  }

  return complex<Float>(re,im);

}

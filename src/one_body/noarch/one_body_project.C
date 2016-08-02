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

complex<Float> OneBody::Project(Float* phi)
{

  Float re = 0.0;
  Float im = 0.0;

  for(int i=0; i<2*vol; i+=2) {
    re += psi[i]*phi[i]+psi[i+1]*phi[i+1];
    im += psi[i]*phi[i+1]-psi[i+1]*phi[i];
  }

  return complex<Float>(re,im);

}

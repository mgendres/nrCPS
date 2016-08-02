#include <iostream>
#include <math.h>
using namespace std;
#include "lattice.h"
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "hamiltonian.h"
#include "dispersion.h"
#include "constants.h"
#include "fourier.h"

void Interaction::Initialize() {

  int vol = GJP.Vol();

  Float* lat = lattice->Get();
  fftComplex* fld = field->Get(0);

  //---- First apply interaction onto lattice
  int j=0;
  for (int i=0; i<vol; ++i) {
    fld[i][0] = lat[j];
    fld[i][1] = lat[j+1];
    j += 2;
  }

//  if ( GJP.APBCQuery() ) { Fourier::phase.Subtract(Fourier::in); }

  field->Backward();

  for (int i=0; i<vol; ++i) {
    fld[i][0] *= interaction[i]; 
    fld[i][1] *= interaction[i]; 
  }

  field->Forward();
//  if ( GJP.APBCQuery() ) { Fourier::phase.Add(Fourier::out); }

  for (int i=0; i<vol; ++i) {
    fld[i][0] += 1.0; 
    fld[i][0] /= vol; 
    fld[i][1] /= vol; 
  }

}

void Potential::Initialize() {}

void Kinetic::Initialize() {}


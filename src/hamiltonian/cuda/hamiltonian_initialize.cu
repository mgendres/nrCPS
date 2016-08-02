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

#include "cuda_utils.h"
#include "cuda_kernels.h"

void Interaction::Initialize() {

  int vol = GJP.Vol();

  Complex_Eq_cFloat<<<BLOCKS,THREADS>>>(field->Get(0), lattice->Get(), vol);
  
//  if ( GJP.APBCQuery() ) { Fourier::phase.Subtract(Fourier::in); }

  field->Backward();

  Complex_TimesEq_Float<<<BLOCKS,THREADS>>>(field->Get(0), dev_interaction, vol);

  field->Forward();

//  if ( GJP.APBCQuery() ) { Fourier::phase.Add(Fourier::out); }

  Complex_PlusEq_One_Div_Vol<<<BLOCKS,THREADS>>>(field->Get(0), vol);

}

void Potential::Initialize() {}

void Kinetic::Initialize() {}

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


void Interactions::Evolve(Fourier* b)
{

  const char* fname = "Interactions::Evolve(Fourier*)";
  VRB.Func(fname);

  ERR.NotImplemented(fname,"GPU implementation is not yet available.");

/*
  for (int j=0; j<2*vol; ++j) { c[j] = 0.0; }
  
  for (int i=0; i<interactions.size(); ++i) {
    for (int j=0; j<2*vol; ++j) { a[j] = b[j]; }
    interactions[i]->Evolve(a);
    for (int j=0; j<2*vol; ++j) { c[j] += a[j]; }
  }

  for (int j=0; j<2*vol; ++j) { b[j] *= 1.0-interactions.size(); }

  for (int j=0; j<2*vol; ++j) { b[j] += c[j]; }
*/
}


void Interaction::Evolve(Fourier* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  const char* fname = "void Interaction3::Evolve(Fourier*)";
  VRB.Func(fname);

  int vol = GJP.Vol();

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"APBCs not permitted."); }

  b->Forward(); //---- Go to position space

//  if ( GJP.APBCQuery() ) { Fourier::phase.Add(Fourier::out); } //---- Add APBC phases...

  for (int n=0; n<b->BatchSize(); ++n) {
    Complex_TimesEq_Complex<<<BLOCKS,THREADS>>>(b->Get(n), field->Get(0), vol );
  }

//  if ( GJP.APBCQuery() ) { Fourier::phase.Subtract(Fourier::in); } //---- Remove APBC phases...

  b->Backward(); //---- Go to momentum space

}

void Potential::Evolve(Fourier* b)
{
  // Apply opererator onto in and return as out (i.e., b -> potential.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  const char* fname = "void Potential::Evolve(Fourier*)";
  VRB.Func(fname);

  int vol = GJP.Vol();

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"APBCs not permitted."); }

  b->Forward(); //---- Go to position space

  for (int n=0; n<b->BatchSize(); ++n) {
    Complex_TimesEq_Float<<<BLOCKS,THREADS>>>(b->Get(n), dev_potential, vol );
  }

  b->Backward(); //---- Go to momentum space

}

void Kinetic::Evolve(Fourier* b)
{
  // Apply opererator onto in and return as out (i.e., b -> kinetic.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  const char* fname = "void Kinetic::Evolve(Fourier*)";
  VRB.Func(fname);

  int vol = GJP.Vol();

  for (int n=0; n<b->BatchSize(); ++n) {
    Complex_DivEq_Float<<<BLOCKS,THREADS>>>(b->Get(n), dev_xi, vol );
  }

}


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

void Interactions::Evolve(Fourier* b)
{

  const char* fname = "Interactions::Evolve(Fourier*)";
  VRB.Func(fname);
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

  const char* fname = "void Interaction::Evolve(Fourier*)";
  VRB.Func(fname);

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"APBCs no allowed."); }

  int vol = GJP.Vol();
  fftComplex* a;
  fftComplex* fld = field->Get(0);

  b->Forward(); //---- Go to position space

//  if ( GJP.APBCQuery() ) { Fourier::phase.Add(Fourier::out); } //---- Add APBC phases

  Float re, im;
  for (int n=0; n<b->BatchSize(); ++n) {
    a = b->Get(n);
    for (int i=0; i<vol; ++i) {
      re = fld[i][0] * a[i][0] - fld[i][1] * a[i][1];
      im = fld[i][0] * a[i][1] + fld[i][1] * a[i][0];
      a[i][0] = re;
      a[i][1] = im;
    }
  }

//  if ( GJP.APBCQuery() ) { Fourier::phase.Subtract(Fourier::in); } //---- Remove APBC phases

  b->Backward(); //---- Go to momentum space

}

void Potential::Evolve(Fourier* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.a ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  const char* fname = "void Potential::Evolve(Fourier*)";
  VRB.Func(fname);

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"APBCs no allowed."); }

  int vol = GJP.Vol();
  fftComplex* a;

  b->Forward(); //---- Go to position space

  for (int n=0; n<b->BatchSize(); ++n) {
    a = b->Get(n);
    for (int i=0; i<vol; ++i) {
      a[i][0] *= potential[i];
      a[i][1] *= potential[i];
    }
  }

  b->Backward(); //---- Go to momentum space

}

void Kinetic::Evolve(Fourier* b)
{
  // Apply opererator onto in and return as out (i.e., b -> interaction.b ).
  // In and out fermion vectors are in the fermion storage order: [x][y][z][comp]

  const char* fname = "void Kinetic::Evolve(Fourier*)";
  VRB.Func(fname);

  int vol = GJP.Vol();
  fftComplex* a;

  //---- Apply D^-1 to source b; b ->  D^-1 b
  for (int n=0; n<b->BatchSize(); ++n) {
    a = b->Get(n);
    for (int i=0; i<vol; ++i) {
      a[i][0] /= xi[i];
      a[i][1] /= xi[i];
    }
  }

}

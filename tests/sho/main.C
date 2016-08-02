#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "constants.h"
#include "enum.h"
#include "arg.h"
#include "lattice.h"
#include "propagator.h"
#include "random.h"
#include "one_body.h"
#include "slater_det.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"

int main()
{

  const char* fname = "int main()";
  cout << "Here we go!\n\n";

  //---- Import simulation parameters

  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");

  DoArg do_arg;
  do_arg.Decode("args/do.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  PropagatorArg propagator_arg;
  propagator_arg.Decode("args/propagator.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  //---- Perform some checks
  if (one_body_arg.source_type!=SOURCE_TYPE_SHO) {
    VRB.Warn(fname,"Setting one_body_arg.source_type to SHO.");
    one_body_arg.source_type=SOURCE_TYPE_SHO;
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);

#define FOURTHORDER
#ifndef FOURTHORDER

  // This is a second order integrator: dt errore are O(dt^3)

  Kinetic kinetic(kinetic_arg,0.5);
  Potential potential(potential_arg,1.0);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

#else

  // This is a fourth order integrator: dt errore are O(dt^5)

  Float a1 =  0.675603595979829;
  Float a2 =  1.35120719195966;
  Float a3 = -0.175603595979829;
  Float a4 = -1.70241438391932;
  Kinetic kinetic_1(kinetic_arg, a1);
  Potential potential_2(potential_arg, a2);
  Kinetic kinetic_3(kinetic_arg, a3);
  Potential potential_4(potential_arg, a4);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic_1);
  hamiltonians.push_back(&potential_2);
  hamiltonians.push_back(&kinetic_3);
  hamiltonians.push_back(&potential_4);
  hamiltonians.push_back(&kinetic_3);
  hamiltonians.push_back(&potential_2);
  hamiltonians.push_back(&kinetic_1);

#endif

  //---- Construct the wavefunctions
  int N = propagator_arg.n_fermions;
  if (N>8) { ERR.General(fname, "propagator_arg.n_fermions must be less than or equal to eight."); }

  OneBody* sources[N];
  char f_name[99];
  for (int n=0; n<N; ++n) {
    sources[n] = new OneBody(one_body_arg);
    sprintf(f_name, "%s.%d", "sources/source", n); 
    sources[n]->Set(f_name);
  }

  //---- Open files for measurments
  //---- Files labeled by particle number
  //---- Each file will contain data for the entire ensemble:
  //---- row = config number, col = time separation
  ofstream file[N];

  for (int n=0;n<N; ++n) {
    sprintf(f_name, "%s.%d", "results/sho", n); 
    file[n].setf(ios_base::scientific);
    file[n].open(f_name,ios_base::app);
  }

  complex<Float> corr;

  //---- Instantiate propagator classes
  Propagator propagator(propagator_arg);

  //--- Set pointers to propagators at time slice t
  Float* props[N];
  for (int n=0; n<N; ++n) {
    props[n] = propagator.Get(n);
  }

  //---- Set propagator sources and inititialize
  for (int n=0; n<N; ++n) {
    propagator.Set(sources[n]->Get(), n);
  }

  //---- Perform contractions and store results for each time slice
  for(int t=0; t<GJP.Tsites(); ++t) {

    VRB.Debug(fname, "t = %d", t);

    //---- Normalize the wave function
//if (t%10==0) {
//  for (int n=0; n<N; ++n) {
//    VRB.Debug(fname, "Normalizing propagator %d of %d\n", n, N);
//    prop[n]->Normalize();
//  }
//} 

    //---- Compute correlation functions
    for (int n=0; n<N; ++n) {
      corr = sources[n]->Project(props[n]);
      file[n] << setprecision(PREC) << corr.real() << " ";
      file[n] << setprecision(PREC) << corr.imag() << " ";
    }

    //---- Compute propagators on the current background configuration
    propagator.Run(hamiltonians);

  }

  for (int n=0; n<N; ++n) {
    file[n] << "\n";
  }

  //---- Delete one-body memory allocation
  for (int n=0; n<N; ++n) {
    delete sources[n];
  }

  //---- Close data files 
  for (int n=0; n<N; ++n) {
    file[n].close();
  }

  return(EXIT_SUCCESS);
}

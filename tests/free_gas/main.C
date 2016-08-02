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
#include "momenta.h"
#include "dispersion.h"
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

  PropagatorArg prop_arg;
  prop_arg.Decode("args/prop.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  MomentaArg momenta_arg;
  momenta_arg.Decode("args/momenta.arg");

  //---- Perform some checks
  if (one_body_arg.source_type != SOURCE_TYPE_MOM) {
    VRB.Warn(fname, "one_body_arg.source_type not supported:");
    VRB.Warn(fname, "\tusing SOURCE_TYPE_MOM instead.");
    one_body_arg.source_type = SOURCE_TYPE_MOM;
  }

  if (momenta_arg.dispersion_type != kinetic_arg.dispersion_type) {
    VRB.Warn(fname, "momenta_arg.dispersion_type and prop_arg.dispersion_type disagree:");
    VRB.Warn(fname, "\tSetting momenta_arg.dispersion_type to prop_arg.dispersion_type.");
    momenta_arg.dispersion_type = kinetic_arg.dispersion_type;
  }

  if (momenta_arg.mass != kinetic_arg.mass) {
    VRB.Warn(fname, "momenta_arg.mass and prop_arg.mass disagree:");
    VRB.Warn(fname, "\tSetting momenta_arg.mass to prop_arg.mass.");
    momenta_arg.mass = kinetic_arg.mass;
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  Kinetic kinetic(kinetic_arg,1.0);
  hamiltonians.push_back(&kinetic);

  //---- Determine lattice momenta below the fermi surface 
  Momenta momenta(momenta_arg);
  int N = momenta.MomentaCount();
  vector<int> momentum;

  //---- Print out some information about momenta
  if (1) {
    Dispersion dispersion(kinetic_arg.dispersion_type, kinetic_arg.mass, CUTOFF_TYPE_HARD);
    cout << "Fermi energy: " << momenta_arg.fermi_energy << "\n";
    cout << "Number of momenta below Fermi surface: " << N << "\n";
    for (int n=0; n<N; ++n) {
      momentum = momenta.GetMomentum(n);
      cout << n << " (" << momenta.OppositeParityIndex(n) << "): ";
      cout << momentum[0] << " " << momentum[1] << " " << momentum[2] << " ";
      cout << log(1+dispersion.Get(momentum[0],momentum[1],momentum[2])) << "\n";
    }
  }

  //---- Some measurement classes for computing Slater determinants of one particle functions
  OneBody* one_bodies[N];
  for (int n=0; n<N; ++n) {
    momentum = momenta.GetMomentum(n);
    one_bodies[n] = new OneBody(one_body_arg);
    one_bodies[n]->Set(momentum[0],momentum[1],momentum[2]);
  }

  //---- Open files for measurments (files labeled by particle number; each file will contain data for the entire ensemble; row = config number, col = time sep.)
  ofstream file[N];

  char f_name[99];
  for (int n=0; n<N; ++n) {
    sprintf(f_name, "%s.%d", "results/slater", n+1); 
    file[n].setf(ios_base::scientific);
    file[n].open(f_name,ios_base::app);
  }

  complex<Float> corr;

  //---- Instantiate propagator classes
  prop_arg.n_fermions = N;
  Propagator propagator(prop_arg);

  //--- Set pointers to propagators at time slice t
  Float* props[N];
  for (int n=0; n<N; ++n) {
    props[n] = propagator.Get(n);
  }

  //---- Set propagator sources and inititialize
  for (int n=0; n<N; ++n) {
    propagator.Set(one_bodies[n]->Get(), n);
  }

  for(int t=0; t<GJP.Tsites(); ++t) {

    VRB.Debug(fname, "t = %d", t);

    //---- Compute Slater determinants of various sub-matrices
    SlaterDet slater_det(props, one_bodies, N);

    for (int n=0; n<N; ++n) {
      corr = slater_det.Run(n+1);
      corr *= conj(corr);
      file[n] << setprecision(PREC) << corr.real() << " " << corr.imag() << "  ";
    }

    //---- Compute propagators on the current background configuration
    propagator.Run(hamiltonians);

  }

  for (int n=0; n<N; ++n) {
    file[n] << "\n";
  }

  //---- Delete one-body memory allocation
  for (int n=0; n<N; ++n) {
    delete one_bodies[n];
  }

  //---- Close data files 
  for (int n=0; n<N; ++n) {
    file[n].close();
  }

  return(EXIT_SUCCESS);
}

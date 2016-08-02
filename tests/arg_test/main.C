#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "enum.h"
#include "arg.h"
#include "global_job_parameter.h"


int main()
{

  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  DoArg do_arg; // Container for lattice parameters
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg);

  LatticeArg lattice_arg; // Container for lattice parameters
  lattice_arg.Decode("args/lattice.arg");

  RandomArg random_arg; // Container for random generator parameters
  random_arg.Decode("args/random.arg");

  PropagatorArg propagator_arg; // Container for propagator parameters
  propagator_arg.Decode("args/propagator.arg");

  TwoBodyArg two_body_arg;  // Container for two-body correlator parameters
  two_body_arg.Decode("args/two_body.arg");

  GeneralizedTwoBodyArg generalized_two_body_arg;  // Container for two-body correlator parameters
  generalized_two_body_arg.Decode("args/generalized_two_body.arg");

  VerboseArg verbose_arg; // Container for verbosity  parameters
  verbose_arg.Decode("args/verbose.arg");

  EvoArg evo_arg; // Container for evolution parameters
  evo_arg.Decode("args/evo.arg");

  MomentaArg momenta_arg; // Container for momenta generation parameters
  momenta_arg.Decode("args/momenta.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  OneBodyArg one_body_arg;  // Container for two-body correlator parameters
  one_body_arg.Decode("args/one_body.arg");

  KineticArg kinetic_arg; // Container for Kinetic parameters
  kinetic_arg.Decode("args/kinetic.arg");

  InteractionArg interaction_arg; // Container for split point interaction
  interaction_arg.Decode("args/interaction.arg");

  //---- Print out current status of Args (used for debugging)
  cout << "\n";
  cout << "Current status of args:" << "\n";
  cout << "\n";

  cout << "do_arg:\n";
  cout << do_arg.t_sites << "\n";
  cout << do_arg.x_sites << "\n";
  cout << do_arg.y_sites << "\n";
  cout << do_arg.z_sites << "\n";
  cout << do_arg.x_boundary_type << "\n";
  cout << do_arg.y_boundary_type << "\n";
  cout << do_arg.z_boundary_type << "\n";
  cout << do_arg.cutoff << "\n";
  cout << "\n";

  cout << "lattice_arg:\n";
  cout << lattice_arg.field_type << "\n";
  cout << "\n";

  cout << "random_arg:\n";
  cout << random_arg.random_type << "\n";
  cout << random_arg.seed_type << "\n";
  cout << random_arg.seed << "\n";
  cout << random_arg.file_stem << "\n";
  cout << "\n";

  cout << "propagator_arg:\n";
  cout << propagator_arg.n_fermions << "\n";
  cout << propagator_arg.file_stem << "\n";
  cout << "\n";

  cout << "two_body_arg:\n";
  cout << two_body_arg.wavefunc_type << "\n";
  cout << two_body_arg.dispersion_type1 << "\n";
  cout << two_body_arg.dispersion_type2 << "\n";
  cout << two_body_arg.mass1 << "\n";
  cout << two_body_arg.mass2 << "\n";
  cout << two_body_arg.lambda << "\n";
  cout << two_body_arg.file_stem << "\n";
  cout << "\n";

  cout << "generalized_two_body_arg:\n";
  cout << generalized_two_body_arg.wavefunc_type << "\n";
  cout << generalized_two_body_arg.distinct_vec.GetX() << " ";
  cout << generalized_two_body_arg.distinct_vec.GetY() << " ";
  cout << generalized_two_body_arg.distinct_vec.GetZ() << " " << endl;
  cout << generalized_two_body_arg.lambda << "\n";
  cout << two_body_arg.file_stem << "\n";
  cout << "\n";

  cout << "verbose_arg:\n";
  cout << verbose_arg.func_level << "\n";
  cout << verbose_arg.warn_level << "\n";
  cout << verbose_arg.result_level << "\n";
  cout << verbose_arg.flow_level << "\n";
  cout << verbose_arg.debug_level << "\n";
  cout << "\n";

  cout << "evo_arg:\n";
  cout << evo_arg.start << "\n";
  cout << evo_arg.unload_period << "\n";
  cout << evo_arg.configurations << "\n";
  cout << "\n";

  
  cout << "momenta_arg:\n";
  cout << momenta_arg.dispersion_type << "\n";
  cout << momenta_arg.fermi_energy << "\n";
  cout << momenta_arg.mass << "\n";
  cout << "\n";

  cout << "kinetic_arg:\n";
  cout << kinetic_arg.dispersion_type << "\n";
  cout << kinetic_arg.mass<< "\n";
  cout << "\n";

  cout << "interaction_arg:\n";
  cout << interaction_arg.mass << endl;
  cout << interaction_arg.interaction_type << "\n";
  for (int n=0; n<interaction_arg.num_couplings; ++n) {
  cout << interaction_arg.couplings[n] << endl;
  }
  cout << interaction_arg.num_couplings << endl;

  cout << "\n";

  cout << "potential_arg:\n";
  cout << potential_arg.potential_form << "\n";
  cout << potential_arg.potential_type << "\n";
  cout << potential_arg.spring_constant1 << "\n";
  cout << potential_arg.spring_constant2 << "\n";
  cout << potential_arg.spring_constant3 << "\n";
  cout << "\n";

  cout << "one_body_arg:\n";
  cout << one_body_arg.source_type << "\n";
  cout << one_body_arg.lambda1 << "\n";
  cout << one_body_arg.lambda2 << "\n";
  cout << one_body_arg.lambda3 << "\n";
  cout << "\n";



  cout << "Boundary Conditions:\n";
  Float bc[3];
  GJP.GetBcShift(bc); 
  cout << "x: " << bc[0] << "\n";
  cout << "y: " << bc[1] << "\n"; 
  cout << "z: " << bc[2] << "\n";
  cout << "\n";

  cout << "GJP.PBCQuery():\n";
  cout << GJP.APBCQuery() << "\n";
  cout << "\n";

  cout << "GJP.Cutoff():\n";
  cout << GJP.Cutoff() << "\n";
  cout << GJP.CutoffSq() << "\n";

  return(EXIT_SUCCESS);
}

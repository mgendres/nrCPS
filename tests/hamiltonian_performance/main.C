#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <time.h>
using namespace std;
#include "constants.h"
#include "enum.h"
#include "arg.h"
#include "lattice.h"
#include "propagator.h"
#include "random.h"
#include "one_body.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "two_body";

  VRB.Flow(fname, "Running production code %s, compiled on %s at %s. ", prog_name, __DATE__, __TIME__ );
  VRB.LibInfo(fname);

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg;
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  LatticeArg lattice2_arg;
  lattice2_arg.Decode("args/lattice2.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  PropagatorArg prop_arg;
  prop_arg.Decode("args/prop.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattice (contains two and three-body fields, initialized when instantiated)
  Lattice lattice2(lattice2_arg, &rng);

  //---- Instantiate two- and three-body hamiltonians, as well as the external potential
  Interaction interaction(&lattice2, interaction_arg, 1.0);
  Potential potential(potential_arg, 0.5);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
//  hamiltonians.push_back(&kinetic);
//  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);

  //---- Some measurement classes for evaluating two and three-body correlators (assumes different species)

  OneBody source(one_body_arg);
  source.Set("args/source.arg");

  //---- Instantiate propagator classes
  Propagator prop(prop_arg);

  //---- Set the propagator sources and inititialize
  for (int n=0; n<prop_arg.n_fermions; ++n) {
    prop.Set(source.Get(), n);
  }

  for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }

  clock_t clock_base = clock();
  Float  clock_elapsed;

  for(int t=0; t<GJP.Tsites(); ++t) {

    //---- Compute propagators on the current background configuration
    prop.Run(hamiltonians);

  }

  clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
  VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

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
#include "hamiltonian.h"
#include "random.h"
#include "one_body.h"
#include "slater_det.h"
#include "verbose.h"
#include "error.h"
#include "momenta.h"
#include "dispersion.h"
#include "global_job_parameter.h"
#include "sysfunc.h"

int main()
{


  //---- Establish communications (needed to pass memory check)
  Comms::Initialize();

  const char* fname = "int main()";
  cout << "Here we go!\n\n";

  //---- Import simulation parameters

  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");

  DoArg do_arg;
  do_arg.Decode("args/do.arg");

  LatticeArg lattice2_arg;
  lattice2_arg.Decode("args/lattice2.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);
  Random rng(random_arg);

  //---- Instantiate the lattice
  Lattice lattice2(lattice2_arg, &rng);

  if (1) {
    //---- Instantiate two- and three-body hamiltonians, as well as the external potential
    Potential potential1(potential_arg, 1.0);
  }

  if (1) {
    //---- Instantiate two- and three-body hamiltonians, as well as the external potential
    Interaction interaction1(&lattice2, interaction_arg, 1.0);
    Potential potential1(potential_arg, 1.0);
    Potential potential2(potential_arg, 1.0);
    Interaction interaction2(&lattice2, interaction_arg, 1.0);
    Potential potential3(potential_arg, 1.0);
  }

  //---- Close communications (needed to pass memory check)
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

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
  VerboseArg verbose_arg("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg("args/do.arg");
  GJP.Initialize(do_arg); 

  LatticeArg lattice_arg("args/lattice.arg");
  EvoArg evo_arg("args/evo.arg");
  RandomArg random_arg("args/random.arg");
  InteractionArg interaction_arg("args/interaction.arg");
  PotentialArg potential_arg("args/potential.arg");
  PropagatorArg prop_arg("args/prop.arg");
  KineticArg kinetic_arg("args/kinetic.arg");
  OneBodyArg one_body_arg("args/one_body.arg");
  TwoBodyArg two_body_arg("args/two_body.arg");

  //---- Perform some checks
  if (evo_arg.unload_period<1) { ERR.General(fname,"evo_arg.unload_period must be greater than zero."); }

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattice (contains two and three-body fields, initialized when instantiated)
  Lattice lattice(lattice_arg, &rng);

  //---- Instantiate two- and three-body hamiltonians, as well as the external potential
  Interaction interaction(&lattice, interaction_arg, 1.0);
//  Potential potential(potential_arg, 0.5);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
//  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
//  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Some measurement classes for evaluating two and three-body correlators (assumes different species)

  OneBody source0(one_body_arg);
  OneBody source1(one_body_arg);
  source0.Set("args/source0.arg");
  source1.Set("args/source1.arg");
  TwoBody two_body(two_body_arg);

  //---- Open files for measurments
  //---- row = config number, col = time separation
  FILE* file;
  char f_name[99];
  sprintf(f_name, "%s", (two_body_arg.file_stem).c_str()); 
  file = Fopen(IO_TYPE_ROOT,f_name,"a");

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr[GJP.Tsites()];
  for (int t=0; t<GJP.Tsites(); ++t) { corr[t]=0.0; }

  //---- Instantiate propagator classes
  if (prop_arg.n_fermions != 2) {ERR.General(fname, "Number of prop_arg.n_fermions must equal two."); }
  Propagator propagator(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  //---- Begin the simulation by looping over configurations
  for (int j=0; j<configurations; ++j) {

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);
    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    //---- Save the current RNG state
    rng.WriteState(random_arg.file_stem.c_str());

    for (int k=0; k<unload_period; ++k) {

      //---- Set the propagator sources and inititialize
      propagator.Set(source0.Get(), 0);
      propagator.Set(source1.Get(), 1);
    
      for(int t=0; t<GJP.Tsites(); ++t) {

//        VRB.Debug(fname, "t = %d", t);
    
        //---- Perform contractions
        corr[t] +=  pow( two_body.Run( propagator.Get(0), propagator.Get(1) ) , 2);

        //---- Compute propagators on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }
    }

    //---- Perform global sums
    Comms::GlobalSum( corr, GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for(int t=0; t<GJP.Tsites(); ++t) {

      //---- Store accumulated results in file 
      Fprintf(file, "%.*e %.*e ", PREC, corr[t].real() / block_size, PREC, corr[t].imag() / block_size);

      //---- Reset accumulated corr
      corr[t] = 0.0;

    }

    Fprintf(file, "\n");

    //--- Perform an occasional flush; this is needed because otherwise if the program quits unexpectedly, all data will be lost
    if (j % FLUSHRATE == 0 ) {
      Fflush(file);
    }

  }

  //---- Close data file
  Fclose(file);

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

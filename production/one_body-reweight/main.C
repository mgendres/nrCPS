#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
#include <time.h>
using namespace std;
#include "constants.h"
#include "enum.h"
#include "arg.h"
#include "lattice.h"
#include "propagator.h"
#include "random.h"
#include "one_body.h"
#include "verbose.h"
#include "error.h"
#include "hamiltonian.h"
#include "global_job_parameter.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "one_body";

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
  PropagatorArg prop_arg("args/prop.arg");
  InteractionArg interaction_arg("args/interaction.arg");
  PotentialArg potential_arg("args/potential.arg");
  KineticArg kinetic_arg("args/kinetic.arg");
  OneBodyArg one_body_arg("args/one_body.arg");

  //---- Perform some checks
  if (evo_arg.unload_period!=1) { ERR.General(fname,"evo_arg.unload_period must be equal to one."); }

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattices for two and three body auxilliary fields
  Lattice lattice(lattice_arg, &rng);

  //---- Instantiate two-body and three-body hamiltonians, as well as the external potential
  Interaction interaction(&lattice, interaction_arg, 1.0);
  Potential potential(potential_arg, 0.5);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Set sources; these sources will be used by Propagator
  OneBody source(one_body_arg);
  source.Set("args/source.arg");
  OneBody sink(one_body_arg);
  sink.Set("args/sink.arg");

  //---- Open files for measurments
  //---- row = config number, col = time separation
  FILE* file;
  char f_name[99];
  sprintf(f_name, "%s", (prop_arg.file_stem).c_str());
  file = Fopen(IO_TYPE_ROOT,f_name,"a");

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr[GJP.Tsites()];
  for (int t=0; t<GJP.Tsites(); ++t) { corr[t]=0.0; }

  //---- Instantiate propagator classs
  if (prop_arg.n_fermions != 1) {ERR.General(fname, "Number of prop_arg.n_fermions must equal one."); }
  Propagator prop(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  Float phi0 = 0.0022;
  Float norm;

  //---- Begin the simulation by looping over configurations
  for (int j=0; j<configurations; ++j) {

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);
    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    //---- Save the current RNG state
    rng.WriteState(random_arg.file_stem.c_str());

    for (int k=0; k<unload_period; ++k) {

      //---- Set propagator sources and inititialize
      prop.Set(source.Get(), 0);

      for(int t=0; t<GJP.Tsites(); ++t) {

//        VRB.Debug(fname, "t = %d", t);

        //---- Perform contractions and store results at each time slice
        corr[t] += sink.Project(prop.Get(0));

        //---- Generate new field configurations
        lattice.Refresh(phi0);
      
        //---- Compute propagator on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        prop.Run(hamiltonians);

        // Normalize the propagators
        norm =   exp( GJP.Vol() * phi0 * ( real( lattice.Mean() ) - 0.5 * phi0) ) ;
        prop.Normalize( norm );
       

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

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

bool fileQ(const char *filename)
{
  ifstream ifile(filename);
  return ifile;
}

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "few_body";

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

  LatticeArg lattice_arg;
  lattice_arg.Decode("args/lattice.arg");

  EvoArg evo_arg;
  evo_arg.Decode("args/evo.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  PropagatorArg prop_arg;
  prop_arg.Decode("args/propagator.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  // Determine the number of two_body.arg files
  int n_two_body=0;
  {
    bool file_exists=1;
    while (file_exists) {
      char f_name[99];
      sprintf(f_name, "args/two_body_%d.arg", n_two_body ); 
      if ( fileQ(f_name) ) {
        n_two_body++;
      } else {
        file_exists=0;
      }
    }
  }

  //---- Open two_body.arg files
  TwoBodyArg two_body_arg[n_two_body];
  for (int b=0; b<n_two_body; ++b) {
    char f_name[99];
    sprintf(f_name, "args/two_body_%d.arg", b); 
    two_body_arg[b].Decode(f_name);
  }

  //---- Perform some checks
  if (evo_arg.unload_period<1) { ERR.General(fname,"evo_arg.unload_period must be greater than zero."); }

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattices
  Lattice lattice(lattice_arg, &rng);

  //---- Instantiate Hamiltonians
  Interaction interaction(&lattice, interaction_arg, 1.0);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to Hamiltonians; these pointers are passed to Propagator.Run(...)
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&interaction);
  hamiltonians.push_back(&kinetic);

  //---- Some measurement classes 
  int N = 3;
  OneBody* one_bodies[N];
  for (int n=0; n<N; ++n) {
    one_bodies[n] = new OneBody(one_body_arg);
  }
  one_bodies[0]->Set(0, 0, 0);                 //  ( 0,  0,  0)
  one_bodies[1]->Set(0, 0, 1);                 //  ( 0,  0,  1)
  one_bodies[2]->Set(0, 0, GJP.Zsites()-1);    //  ( 0,  0, -1)

  //----
  TwoBody* two_bodies[n_two_body];
  for (int b=0; b<n_two_body; ++b) {
    two_bodies[b] = new TwoBody(two_body_arg[b]);
    two_bodies[b]->Deform("args/two_body_def.arg");
  }

  //---- Open files for measurments
  //---- row = config number, col = time separation
  FILE* file[n_two_body];
  for (int b=0; b<n_two_body; ++b) {
    char f_name[99];
    sprintf(f_name, "%s", (two_body_arg[b].file_stem).c_str()); 
    file[b] = Fopen(IO_TYPE_ROOT,f_name,"a");
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;

  complex<Float> corr[n_two_body][GJP.Tsites()];
  for (int b=0; b<n_two_body; ++b) {
    for (int t=0; t<GJP.Tsites(); ++t) { corr[b][t] = 0.0; }
  }

  //---- Instantiate propagator classes
  if (prop_arg.n_fermions != N) {
    VRB.Warn(fname, "Setting prop_arg.n_fermions equal to %d.", N);
    prop_arg.n_fermions = N;
  }
  Propagator propagator(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  complex<Float> two00, two12, two10, two02;

  //---- Begin the simulation by looping over configurations
  for (int j=0; j<configurations; ++j) {

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);
    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    //---- Save the current RNG state
    rng.WriteState(random_arg.file_stem.c_str());

    for (int k=0; k<unload_period; ++k) {

      //---- Set the propagator sources and inititialize
      for (int n=0; n<N; ++n) {
        propagator.Set(one_bodies[n]->Get(), n);
      }

      for(int t=0; t<GJP.Tsites(); ++t) {

        for(int b=0; b<n_two_body; ++b) {
          two00 = two_bodies[b]->Run(propagator.Get(0), propagator.Get(0));
          two12 = two_bodies[b]->Run(propagator.Get(1), propagator.Get(2));
          two10 = two_bodies[b]->Run(propagator.Get(1), propagator.Get(0));
          two02 = two_bodies[b]->Run(propagator.Get(0), propagator.Get(2));

          //----
          //---- Perform contractions for three fermions
          //----

          // u: [0] = ( 0, 0, 0) and [1] = ( 0, 0, 1)
          // d: [0] = ( 0, 0, 0) and [2] = ( 0, 0,-1)
          corr[b][t] += two00*two12 - two10*two02;
  
        }

        //---- Compute propagators on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }
    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr[0][0]), n_two_body * GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for (int b=0; b<n_two_body; ++b ) {
      for(int t=0; t<GJP.Tsites(); ++t) {
        //---- Store accumulated results in file 
        Fprintf(file[b], "%.*e %.*e ", PREC, corr[b][t].real() / block_size, PREC, corr[b][t].imag() / block_size);
        //---- Reset accumulated corr
        corr[b][t] = 0.0;
      }

      Fprintf(file[b], "\n");

      //---- Perform an occasional flush; 
      //---- this is needed because otherwise if the program quits unexpectedly, all data will be lost
      if (j % FLUSHRATE == 0 ) {
        Fflush(file[b]);
      }
    }
  }

  //---- Close data file
  //---- Delete wave function body memory allocation
  for (int n=0; n<N; ++n) {
    delete one_bodies[n];
  }
  for (int b=0; b<n_two_body; ++b) {
    delete two_bodies[b];
    Fclose(file[b]);
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

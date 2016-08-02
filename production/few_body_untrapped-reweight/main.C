#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <time.h>
#include <math.h>
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
  const char* prog_name = "few_body";

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
  PropagatorArg prop_arg("args/propagator.arg");
  KineticArg kinetic_arg("args/kinetic.arg");
  OneBodyArg one_body_arg("args/one_body.arg");
  TwoBodyArg two_body_arg("args/two_body.arg");

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
  TwoBody two_body(two_body_arg);
  two_body.Deform("args/two_body_def.arg");

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
  for (int t=0; t<GJP.Tsites(); ++t) { corr[t] = 0.0; }

  //---- Instantiate propagator classes
  if (prop_arg.n_fermions != N) {
    VRB.Warn(fname, "Setting prop_arg.n_fermions equal to %d.", N);
    prop_arg.n_fermions = N;
  }

  Propagator propagator(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  complex<Float> two00, two12, two10, two02;;

  Float phi0;

//Float x0= 0.0014;
//Float y0= 0.0000;

//Float x0= 0.0;
//Float y0= 0.0081;

  Float x0= 0.0000;
  Float y0= 0.0081;

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

      phi0=0.0081;
      lattice.Refresh(phi0);
  
      for(int t=0; t<GJP.Tsites(); ++t) {

        two00 = two_body.Run(propagator.Get(0), propagator.Get(0));
        two12 = two_body.Run(propagator.Get(1), propagator.Get(2));
        two10 = two_body.Run(propagator.Get(1), propagator.Get(0));
        two02 = two_body.Run(propagator.Get(0), propagator.Get(2));

        corr[t] += two00*two12 - two10*two02;

        //---- Compute propagators on the current background configuration
        propagator.Normalize(  exp( 0.25 * GJP.Vol() * phi0 * ( real( lattice.Mean() ) - 0.5 * phi0) )  );
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        //---- Generate a new field configuration


        phi0 = x0*log( real( corr[t]) ) + y0;
//        cout << "t = "<< t << ", :" << phi0 << endl;
        lattice.Refresh(phi0);

      }
    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr[0]), GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for(int t=0; t<GJP.Tsites(); ++t) {
      //---- Store accumulated results in file 
      Fprintf(file, "%.*e %.*e ", PREC, corr[t].real() / block_size, PREC, corr[t].imag() / block_size);
      //---- Reset accumulated corr
      corr[t] = 0.0;
    }

    Fprintf(file, "\n");

    //---- Perform an occasional flush; 
    //---- this is needed because otherwise if the program quits unexpectedly, all data will be lost
    if (j % FLUSHRATE == 0 ) {
      Fflush(file);
    }
  }

  //---- Close data file
  Fclose(file);

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

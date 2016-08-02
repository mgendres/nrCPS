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
#include "utils.h"

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

  Float phi0[ 2*GJP.Vol() ];
  Float norm;

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

        // Compute correlator
        two00 = two_body.Run(propagator.Get(0), propagator.Get(0));
        two12 = two_body.Run(propagator.Get(1), propagator.Get(2));
        two10 = two_body.Run(propagator.Get(1), propagator.Get(0));
        two02 = two_body.Run(propagator.Get(0), propagator.Get(2));
        corr[t] += two00*two12 - two10*two02;

        // Generate a new phi0
        if (0) {
          for (int i=0; i<GJP.Vol(); ++i) { phi0[2*i] = 0.0000; phi0[2*i+1] = 0.0000; }
        }
        if (0) {
          for (int i=0; i<GJP.Vol(); ++i) { phi0[2*i] = 0.0081; phi0[2*i+1] = 0.0000; }
        }
        if (0) {
          complex<Float> c_old(0.0,0.0);
          complex<Float> c_new(0.0,0.0);
          Propagator propagator2(prop_arg);
          for (int i=0; i< 2*GJP.Vol(); ++i) { phi0[i] = 0.00; }
          for (int u=0; u<10; ++u) {
            for (int m=0; m<N; ++m) (propagator2.Set(propagator.Get(m), m));
            lattice.Refresh(0.00);
            for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
            propagator2.Run(hamiltonians);
            two00 = two_body.Run(propagator2.Get(0), propagator2.Get(0));
            two12 = two_body.Run(propagator2.Get(1), propagator2.Get(2));
            two10 = two_body.Run(propagator2.Get(1), propagator2.Get(0));
            two02 = two_body.Run(propagator2.Get(0), propagator2.Get(2));
            c_new = two00*two12 - two10*two02;
            c_new /= exp(- 0.5 * InnerProduct( &(lattice[0]), &(lattice[0]), 2*GJP.Vol() ));
            cout << abs(c_old) << " " << abs(c_new) << endl;
            if (abs(c_new)>abs(c_old)) {
              for (int ii=0; ii< 2*GJP.Vol(); ++ii) { phi0[ii] = lattice[ii]; }
              c_old=c_new;
            } 
          }
        }

        //---- Generate a new field configuration
        lattice.Refresh(phi0);

        //---- Compute propagators on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        // Normalize the propagators
        norm =  exp( 0.25 * ( InnerProduct( &(lattice[0]), phi0, 2*GJP.Vol() )  - 0.5 * InnerProduct(phi0, phi0, 2*GJP.Vol() ) ) );
        propagator.Normalize( norm );

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

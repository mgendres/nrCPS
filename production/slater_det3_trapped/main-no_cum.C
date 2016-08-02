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
#include "slater_det2.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "slater_det2_trapped";

  VRB.Flow(fname, "Running production code %s, compiled on %s at %s. ", prog_name, __DATE__, __TIME__ );
  VRB.LibInfo(fname);

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Import simulation parameters

  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");

  DoArg do_arg;
  do_arg.Decode("args/do.arg");

  LatticeArg lattice_arg;
  lattice_arg.Decode("args/lattice.arg");

  EvoArg evo_arg;
  evo_arg.Decode("args/evo.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  PropagatorArg prop_arg;
  prop_arg.Decode("args/propagator.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  //---- Perform some checks
  if (evo_arg.unload_period<1) {
    ERR.General(fname,"evo_arg.unload_period must be greater than zero.");
  }
  if (one_body_arg.source_type!=SOURCE_TYPE_SHO) {
    VRB.Warn(fname,"Setting one_body_arg.source_type to SOURCE_TYPE_SHO.");
    one_body_arg.source_type=SOURCE_TYPE_SHO;
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);
  Random rng(random_arg);

  //---- Instantiate the lattices
  Lattice lattice(lattice_arg, &rng);

//---- Uncomment the following line to use O(dt^2) imroved time discretization
//---- otherwise the O(dt) discretization will be used
#define IMPROVED 
#if defined(IMPROVED)

  VRB.Flow(fname,"Using O(dt^2) improved discretization");

  //---- Instantiate hamiltonians and external potential
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

#else

  VRB.Flow(fname,"Using O(dt) discretization");

  //---- Instantiate hamiltonians and external potential
  Interaction interaction(&lattice, interaction_arg, 1.0);
  Potential potential(potential_arg, 1.0);
  Kinetic kinetic(kinetic_arg, 0.5);

  vector<Hamiltonian*> interaction_list;
  interaction_list.push_back(&interaction);
  interaction_list.push_back(&potential);
  Interactions interactions(interaction_list);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&interactions);
  hamiltonians.push_back(&kinetic);

#endif


  //---- Construct the wavefunctions
  int N = prop_arg.n_fermions;
  if (N>35) { ERR.General(fname, "propagator_arg.n_fermions must be less than or equal to 35."); }
  
  OneBody* one_bodies[N];
  char f_name[99];

  for (int n=0; n<N; ++n) {
    one_bodies[n] = new OneBody(one_body_arg);
    sprintf(f_name, "%s.%d", "sources/source", n); 
    one_bodies[n]->Set(f_name);
  }

  TwoBody two_body(two_body_arg);

  //---- Open files for measurments
  //---- Files labeled by particle number
  //---- Each file will contain data for the entire ensemble:
  //---- row = config number, col = time separation
  FILE* file[2*N];
  for (int n=0; n<2*N; ++n) {
    //---- This is for the single particle propagators
    sprintf(f_name, "%s.%d", "results/slater", n+1); 
    file[n] = Fopen(IO_TYPE_ROOT,f_name,"a");
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr[2*N][GJP.Tsites()];
  for (int n=0; n<2*N; ++n)
  for (int t=0; t<GJP.Tsites(); ++t) {
      corr[n][t]=0.0;
  }

  //---- Instantiate propagator classes
  Propagator propagator(prop_arg);

  //--- Set pointers to propagators at time slice t
  Float* props[N];
  for (int n=0; n<N; ++n) {
    props[n] = propagator.Get(n);
  }

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

      //---- Set propagator sources and inititialize
      for (int n=0; n<N; ++n) {
        propagator.Set(one_bodies[n]->Get(), n);
      }

      for(int t=0; t<GJP.Tsites(); ++t) {

        //VRB.Debug(fname, "t = %d", t);

        //---- Compute Slater determinants of various sub-matrices
        SlaterDet2 slater_det(props, props, one_bodies, NULL,  &two_body, N);

        for (int n=0; n<N; ++n) {
          //---- Case 1: n_up = n+1 , n_down = n+1
          corr[2*n+1][t] += slater_det.Run(n+1, -1, -1);
          //---- Case 2: n_up = n+1, n_down = n
          corr[2*n][t] += slater_det.Run(n+1, -1, n); 
        }

        //---- Compute propagators on the current background configuration
        for (int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }

    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr[0][0]), 2*N * GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for (int n=0; n<2*N; ++n) {
      for(int t=0; t<GJP.Tsites(); ++t) {

        //---- Store accumulated results in file 
        Fprintf(file[n], "%.*e %.*e ", PREC, corr[n][t].real() / block_size, PREC, corr[n][t].imag() / block_size);

        //---- Reset accumulated corr
        corr[n][t]=0.0;
      }

      Fprintf(file[n], "\n");

      //--- Perform an occasional flush; this is needed because otherwise if the program quits unexpectedly, all data will be lost
      if (j % FLUSHRATE == 0 ) {
        Fflush(file[n]);
      }

    }

  }

  //---- Delete one-body memory allocation
  for (int n=0; n<N; ++n) {
    delete one_bodies[n];
  }

  //---- Close data files 
  for (int n=0; n<2*N; ++n) {
    Fclose(file[n]);
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

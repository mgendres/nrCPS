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

  PotentialArg potential_arg;
  potential_arg.Decode("args/potential.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  //---- Perform some checks
  if (one_body_arg.source_type!=SOURCE_TYPE_MOM) { ERR.General(fname,"SourceType not permitted."); }
  if (evo_arg.unload_period<1) { ERR.General(fname,"evo_arg.unload_period must be greater than zero."); }

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattices
  Lattice lattice(lattice_arg, &rng);

  //---- Instantiate Hamiltonians
  Interaction interaction(&lattice, interaction_arg, 1.0);
  Kinetic kinetic(kinetic_arg, 0.5);
  Potential potential(potential_arg, 0.5);

  //---- Create pointers to Hamiltonians; these pointers are passed to Propagator.Run(...)
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Some measurement classes 
  int n_reps = 10;
  OctahedralRep rep[n_reps];
  Vector3 distinct_vecs[n_reps];
  rep[0]  = OCTAHEDRAL_REP_A1g;
  distinct_vecs[0] = Vector3(0,0,0);
  rep[1]  = OCTAHEDRAL_REP_Eg;
  distinct_vecs[1] = Vector3(0,0,1);
  rep[2]  = OCTAHEDRAL_REP_T1u;
  distinct_vecs[2] = Vector3(0,0,1);
  rep[3]  = OCTAHEDRAL_REP_T2g;
  distinct_vecs[3] = Vector3(0,1,1);
  rep[4]  = OCTAHEDRAL_REP_T2u;
  distinct_vecs[4] = Vector3(0,1,1);
  rep[5]  = OCTAHEDRAL_REP_A2u;
  distinct_vecs[5] = Vector3(1,1,1);
  rep[6]  = OCTAHEDRAL_REP_A2g;
  distinct_vecs[6] = Vector3(0,1,2);
  rep[7]  = OCTAHEDRAL_REP_T1g;
  distinct_vecs[7] = Vector3(0,1,2);
  rep[8]  = OCTAHEDRAL_REP_Eu;
  distinct_vecs[8] = Vector3(1,1,2);
  rep[9]  = OCTAHEDRAL_REP_A1u;
  distinct_vecs[9] = Vector3(1,2,3);

  // Have two options: either R = A1g for the single fermion correlator
  OneBody* one_bodies_R[n_reps];
  for (int j=0; j<n_reps; ++j) {
    one_bodies_R[j] = new OneBody(one_body_arg);
    one_bodies_R[j]->Set(-1*distinct_vecs[j], rep[j]); 
  }

  int N = n_reps; // Number of reps simulated
  //int N = 6; // Number of reps simulated

  //---- Open files for measurments
  //---- row = config number, col = time separation
  FILE* file[N];
  for (int b=0; b<N; ++b) {
    char f_name[99];
    sprintf(f_name, "results/corr.%d", b); 
    file[b] = Fopen(IO_TYPE_ROOT,f_name,"a");
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;

  complex<Float> corr[N][GJP.Tsites()];
  for (int b=0; b<N; ++b) {
    for (int t=0; t<GJP.Tsites(); ++t) { corr[b][t] = 0.0; }
  }

  //---- Instantiate propagator classes
  if (prop_arg.n_fermions != N) {
    VRB.Warn(fname, "Setting prop_arg.n_fermions equal to %d.", N);
    prop_arg.n_fermions = N;
  }
  Propagator propagator_A1g(prop_arg);
  Propagator propagator_R(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  complex<Float> oneA, oneB;

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
        propagator_R.Set(one_bodies_R[n]->Get(), n);
      }

      for(int t=0; t<GJP.Tsites(); ++t) {

        for (int b=0; b<N; ++b) {

          // u: [0] = ( 0, 0, 0) and [1] = ( 0, 0, 1)
          // d: [0] = ( 0, 0, 0) and [2] = ( 0, 0,-1)

          corr[b][t] += one_bodies_R[b]->Project( propagator_R.Get(b) );

        }

        //---- Compute propagators on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator_R.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }
    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr[0][0]), N * GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for (int b=0; b<N; ++b ) {
      for(int t=0; t<GJP.Tsites(); ++t) {
        //---- Store accumulated results in file 
        Fprintf(file[b], "%.*e %.*e ", PREC, corr[b][t].real() / block_size, PREC, corr[b][t].imag() / block_size);
        //---- Reset accumulated corr
        corr[b][t] = 0.0;
      }
      Fprintf(file[b], "\n");

      //---- Perform an occasional flush; 
      //---- this is needed because otherwise if the program quits unexpectedly, all data will be lost
      if (j % FLUSHRATE == 0 ) { Fflush(file[b]); }
    }
  }

  //---- Close data file
  //---- Delete wave function body memory allocation
  for (int b=0; b<N; ++b) {
    delete one_bodies_R[b];
    Fclose(file[b]);
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

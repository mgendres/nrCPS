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

  //---- Open two_body.arg files
  GeneralizedTwoBodyArg two_body_arg("args/two_body.arg");

  //---- Perform some checks
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
//  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
//  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Some measurement classes 
  OctahedralRep rep[10];
  int rep_vec_index[10];
  rep[0] = OCTAHEDRAL_REP_A1g; rep_vec_index[0] = 0;
  rep[1] = OCTAHEDRAL_REP_Eg ; rep_vec_index[1] = 1;
  rep[2] = OCTAHEDRAL_REP_T1u; rep_vec_index[2] = 1;
  rep[3] = OCTAHEDRAL_REP_T2g; rep_vec_index[3] = 2;
  rep[4] = OCTAHEDRAL_REP_T2u; rep_vec_index[4] = 2;
  rep[5] = OCTAHEDRAL_REP_A2u; rep_vec_index[5] = 3;
  rep[6] = OCTAHEDRAL_REP_A2g; rep_vec_index[6] = 4;
  rep[7] = OCTAHEDRAL_REP_T1g; rep_vec_index[7] = 4;
  rep[8] = OCTAHEDRAL_REP_Eu ; rep_vec_index[8] = 5;
  rep[9] = OCTAHEDRAL_REP_A1u; rep_vec_index[9] = 6;

  Vector3 distinct_vecs[7];
  distinct_vecs[0] = Vector3(0,0,0);
  distinct_vecs[1] = Vector3(0,0,1);
  distinct_vecs[2] = Vector3(0,1,1);
  distinct_vecs[3] = Vector3(1,1,1);
  distinct_vecs[4] = Vector3(0,1,2);
  distinct_vecs[5] = Vector3(1,1,2);
  distinct_vecs[6] = Vector3(1,2,3);

  // Define sources and sinks
  int n_reps = 3;
  int n_vecs = 2;

  // Sources
  Vector3 shift(1,1,1);
  OneBody* sourceA[n_vecs];
  OneBody* sourceB[n_vecs];
  for (int b=0; b<n_vecs; ++b) {
    sourceB[b] = new OneBody(one_body_arg);
    sourceB[b]->Set( distinct_vecs[b] );
    sourceA[b] = new OneBody(one_body_arg);
    sourceA[b]->Set( shift - distinct_vecs[b] );
  }
  OneBody sourceC(one_body_arg);
  sourceC.Set(-1*shift);

  // Sinks
  GeneralizedTwoBody* two_bodies[n_vecs];
  OneBody* one_bodies[n_vecs];
  for (int b=0; b<n_vecs; ++b) {
    two_body_arg.distinct_vec = distinct_vecs[b];
    two_bodies[b] = new GeneralizedTwoBody(two_body_arg);
    one_bodies[b] = new OneBody(one_body_arg);
    one_bodies[b]->Set(  distinct_vecs[ rep_vec_index[b] ] , OCTAHEDRAL_REP_A1g);
  }


  //---- Open files for measurments
  //---- row = config number, col = time separation
  FILE* file3[n_reps];
  FILE* file4[n_reps];
  for (int r=0; r<n_reps; ++r) {
    char f_name[99];
    sprintf(f_name, "%s_%d.3", (two_body_arg.file_stem).c_str(), r); 
    file3[r] = Fopen(IO_TYPE_ROOT,f_name,"a");
    sprintf(f_name, "%s_%d.4", (two_body_arg.file_stem).c_str(), r); 
    file4[r] = Fopen(IO_TYPE_ROOT,f_name,"a");
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;

  complex<Float> corr3[n_reps][GJP.Tsites()];
  complex<Float> corr4[n_reps][GJP.Tsites()];
  for (int r=0; r<n_reps; ++r) {
    for (int t=0; t<GJP.Tsites(); ++t) {
      corr3[r][t] = 0.0;
      corr4[r][t] = 0.0;
    }
  }

  //---- Instantiate propagator classes
  prop_arg.n_fermions = n_vecs;
  Propagator propagatorA(prop_arg);
  Propagator propagatorB(prop_arg);
  prop_arg.n_fermions = 1;
  Propagator propagatorC(prop_arg);

  clock_t clock_base = clock();
  Float  clock_elapsed;

  // A, B correspond to opposing momenta
  complex<Float> one_A1g_A, one_A1g_B, one_A1g_C;
  complex<Float> two_A1g_0B, two_A1g_AB, two_A1g_0C, two_A1g_AC;
  complex<Float> two_R_0B, two_R_AB, two_R_0C, two_R_AC;


  //---- Begin the simulation by looping over configurations
  for (int j=0; j<configurations; ++j) {

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);
    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    //---- Save the current RNG state
    rng.WriteState(random_arg.file_stem.c_str());

    for (int k=0; k<unload_period; ++k) {

      //---- Set the propagator sources and inititialize
      for (int b=0; b<n_vecs; ++b) {
        propagatorA.Set(sourceA[b]->Get(), b);
        propagatorB.Set(sourceB[b]->Get(), b);
      }
      propagatorC.Set(sourceC.Get(), 0);


      for(int t=0; t<GJP.Tsites(); ++t) {

        for (int r=0; r<n_reps; ++r) {

          int b = rep_vec_index[r]; // For the sink

          Float* prop0 = propagatorB.Get(0);
          Float* propA = propagatorA.Get(b);
          Float* propB = propagatorB.Get(b);
          Float* propC = propagatorC.Get(0);

          two_A1g_0C = two_bodies[b]->Run(prop0, propC, OCTAHEDRAL_REP_A1g);
          two_A1g_AC = two_bodies[b]->Run(propA, propC, OCTAHEDRAL_REP_A1g);
          two_A1g_0B = two_bodies[b]->Run(prop0, propB, OCTAHEDRAL_REP_A1g);
          two_A1g_AB = two_bodies[b]->Run(propA, propB, OCTAHEDRAL_REP_A1g);

          two_R_0C = two_bodies[b]->Run(prop0, propC, rep[r]);
          two_R_AC = two_bodies[b]->Run(propA, propC, rep[r]);
          two_R_0B = two_bodies[b]->Run(prop0, propB, rep[r]);
          two_R_AB = two_bodies[b]->Run(propA, propB, rep[r]);

	  //one_A1g_A = one_bodies[b]->Project( propA );
	  one_A1g_B = one_bodies[b]->Project( propB );
          one_A1g_C = one_bodies[b]->Project( propC );

          corr3[r][t] += two_R_AB * one_A1g_C - two_R_AC*one_A1g_B;
          corr4[r][t] += two_R_0C * two_A1g_AB - two_R_AC * two_A1g_0B - two_R_0B * two_A1g_AC + two_R_AB * two_A1g_0C;

        }

        //---- Compute propagators on the current background configuration
        for (unsigned int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagatorA.Run(hamiltonians);
        propagatorB.Run(hamiltonians);
        propagatorC.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }
    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr3[0][0]), n_reps * GJP.Tsites() );
    Comms::GlobalSum( &(corr4[0][0]), n_reps * GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for (int r=0; r<n_reps; ++r ) {
      for(int t=0; t<GJP.Tsites(); ++t) {
        //---- Store accumulated results in file 
        Fprintf(file3[r], "%.*e %.*e ", PREC, corr3[r][t].real() / block_size, PREC, corr3[r][t].imag() / block_size);
        Fprintf(file4[r], "%.*e %.*e ", PREC, corr4[r][t].real() / block_size, PREC, corr4[r][t].imag() / block_size);
        //---- Reset accumulated corr
        corr3[r][t] = 0.0;
        corr4[r][t] = 0.0;
      }
      Fprintf(file3[r], "\n");
      Fprintf(file4[r], "\n");

      //---- Perform an occasional flush; 
      //---- this is needed because otherwise if the program quits unexpectedly, all data will be lost
      if (j % FLUSHRATE == 0 ) {
        Fflush(file3[r]);
        Fflush(file4[r]);
      }
    }
  }

  //---- Close data file
  //---- Delete wave function body memory allocation
  for (int r=0; r<n_reps; ++r) {
    Fclose(file3[r]);
    Fclose(file4[r]);
  }

  for (int b=0; b<n_vecs; ++b) {
    delete two_bodies[b];
    delete one_bodies[b];
    delete sourceA[b];
    delete sourceB[b];
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

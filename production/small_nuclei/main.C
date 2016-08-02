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
  const char* prog_name = "small_nuclei";

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

  LatticeArg lattice2_arg;
  lattice2_arg.Decode("args/lattice2.arg");

  EvoArg evo_arg;
  evo_arg.Decode("args/evo.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  InteractionArg interaction_arg;
  interaction_arg.Decode("args/interaction.arg");

  InteractionArg interaction3_arg;
  interaction3_arg.Decode("args/interaction3.arg");

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

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);
  Random rng(random_arg);

  if ( GJP.APBCQuery() ) {
    ERR.General(fname,"Sources are not set up for APBCs; use PBCs instead.");
  }

  //---- Instantiate the lattices
  Lattice lattice(lattice2_arg, &rng); // This is for tau0, z2 interaction 
  Lattice lattice3(lattice2_arg, &rng); // This is for tau3, z2 interaction

  //---- Instantiate hamiltonians
  Kinetic kinetic(kinetic_arg, 0.5);
  Interaction3 interaction(&lattice, interaction_arg, kinetic_arg, kinetic_arg, 1.0);
  Interaction3 interaction3_plus(&lattice3, interaction3_arg, kinetic_arg, kinetic_arg, 1.0);
  Interaction3 interaction3_minus(&lattice3, interaction3_arg, kinetic_arg, kinetic_arg, -1.0);

  //---- Create up fermion interactions
  vector<Hamiltonian*> up_interaction_list;
  up_interaction_list.push_back(&interaction);
  up_interaction_list.push_back(&interaction3_plus);
  Interactions up_interactions(up_interaction_list);

  //---- Create down fermion interactions
  vector<Hamiltonian*> down_interaction_list; 
  down_interaction_list.push_back(&interaction);
  down_interaction_list.push_back(&interaction3_minus);
  Interactions down_interactions(down_interaction_list);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> up_hamiltonians;
  vector<Hamiltonian*> down_hamiltonians;

  up_hamiltonians.push_back(&kinetic);
  up_hamiltonians.push_back(&up_interactions);
  up_hamiltonians.push_back(&kinetic);

  down_hamiltonians.push_back(&kinetic);
  down_hamiltonians.push_back(&down_interactions);
  down_hamiltonians.push_back(&kinetic);

  OneBody one_body(one_body_arg);
  one_body.Set("args/source.arg");

  TwoBody two_body(two_body_arg);
  two_body.Deform("args/two_body_def.arg");

  int N = 8; // This should never be changed

  //---- Open files for measurments
  //---- Files labeled by particle number
  //---- Each file will contain data for the entire ensemble:
  //---- row = config number, col = time separation

  FILE* file[N];
  string filename[N];

  filename[0] = "1S0_1body";
  filename[1] = "3S1_1body";
  filename[2] = "Triton_1body";
  filename[3] = "Alpha_1body";
  filename[4] = "1S0_2body";
  filename[5] = "3S1_2body";
  filename[6] = "Triton_2body";
  filename[7] = "Alpha_2body";

  char f_name[99];
  for (int n=0; n<N; ++n) {
    sprintf(f_name,"results/%s", filename[n].c_str()); 
    file[n] = Fopen(IO_TYPE_ROOT,f_name,"a");
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;

  complex<Float> corr[N][GJP.Tsites()];

  for (int n=0; n<N; ++n)
  for (int t=0; t<GJP.Tsites(); ++t) {
      corr[n][t]=0.0; 
  }

  complex<Float> u_projection;
  complex<Float> d_projection;
  complex<Float> uu_projection;
  complex<Float> ud_projection;
  complex<Float> dd_projection;

  //---- Instantiate propagator classes
  Propagator uprop(prop_arg);
  Propagator dprop(prop_arg);
  Float* uprops = uprop.Get();
  Float* dprops = dprop.Get();

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
      uprop.SetSource(one_body.Get());
      dprop.SetSource(one_body.Get());

      for(int t=0; t<GJP.Tsites(); ++t) {

        //VRB.Debug(fname, "t = %d", t);

        u_projection = one_body.Project(uprops);
        d_projection = one_body.Project(dprops);
        uu_projection = two_body.Run(uprops, uprops);
        ud_projection = two_body.Run(uprops, dprops);
        dd_projection = two_body.Run(dprops, dprops);

        corr[0][t] += u_projection * d_projection;
        corr[1][t] += u_projection * u_projection;
        corr[2][t] += u_projection * u_projection * d_projection;
        corr[3][t] += u_projection * u_projection * d_projection * d_projection;

        corr[4][t] += ud_projection;
        corr[5][t] += uu_projection;
        corr[6][t] += uu_projection * d_projection;
        corr[7][t] += uu_projection * dd_projection;

        //---- Compute propagators on the current background configuration
        kinetic.Initialize();
        interaction.Initialize();
        interaction3_plus.Initialize();
        interaction3_minus.Initialize();
        uprop.Run(up_hamiltonians);
        dprop.Run(down_hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();
        lattice3.Refresh();

      }

    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr[0][0]), N * GJP.Tsites() );

    //---- Write results to files
    Float block_size = Comms::Size() * unload_period;
    for (int n=0; n<N; ++n) {
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

  //---- Close data files 
  //---- Delete data files

  for (int n=0; n<N; ++n) {
    Fclose(file[n]);
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

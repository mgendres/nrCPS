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
#include "slater_det.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "sysfunc.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "slater_mat2";

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
  if (Comms::Size()!=1) {
    ERR.General(fname,"Number of nodes should equal one.");
  }
  if (evo_arg.unload_period!=1) {
    ERR.General(fname,"evo_arg.unload_period should equal one.");
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);
  Random rng(random_arg);

  //---- Instantiate the lattices
  Lattice lattice2(lattice2_arg, &rng);

  //---- Instantiate two- and three-body hamiltonians, as well as the external potential
  Interaction3 interaction(&lattice2, interaction_arg, kinetic_arg, kinetic_arg, 1.0);
  Potential potential(potential_arg, 0.5);
  Kinetic kinetic(kinetic_arg, 0.5);

  //---- Create pointers to hamiltonians; these pointers are passed to Propagator 
  vector<Hamiltonian*> hamiltonians;
  hamiltonians.push_back(&kinetic);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&interaction);
  hamiltonians.push_back(&potential);
  hamiltonians.push_back(&kinetic);

  //---- Construct the wavefunctions
  int N = 7;
  OneBody* sources[N];
  char f_name[99];
  for (int n=0; n<N; ++n) {
    sources[n] = new OneBody(one_body_arg);
    sprintf(f_name, "%s.%d", "sho_sources/source", n); 
    sources[n]->Set(f_name);
  }

  TwoBody two_body(two_body_arg);

  //---- Open files for measurments
  //---- Files labeled by particle number
  //---- Each file will contain data for the entire ensemble:
  //---- row = config number, col = time separation
  ofstream file[N][N];
  ofstream file2[N][N];

  for (int m=0; m<N; ++m) {
    for (int n=0;n<N; ++n) {
      //---- This is for the single particle propagators
      sprintf(f_name, "%s.%d.%d", "results/slater", m, n); 
      file[m][n].setf(ios_base::scientific);
      file[m][n].open(f_name,ios_base::app);
      //---- This is for multi-particle propagators
      sprintf(f_name, "%s.%d.%d", "results/slater2", m, n); 
      file2[m][n].setf(ios_base::scientific);
      file2[m][n].open(f_name,ios_base::app);
    }
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr[N][N][GJP.Tsites()];
  complex<Float> corr2[N][N][GJP.Tsites()];
  for (int m=0; m<N; ++m)
  for (int n=0; n<N; ++n)
  for (int t=0; t<GJP.Tsites(); ++t) {
      corr[m][n][t]=0.0;
      corr2[m][n][t]=0.0;
  }

  //---- Instantiate propagator classes
  Propagator* prop[N];
  for (int n=0; n<N; ++n) {
    prop[n] = new Propagator(prop_arg);
  }

  //--- Set pointers to propagators at time slice t
  Float* props[N];
  for (int n=0; n<N; ++n) {
    props[n] = prop[n]->Get();
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
        prop[n]->SetSource(sources[n]->Get());
      }

      for(int t=0; t<GJP.Tsites(); ++t) {

        //VRB.Debug(fname, "t = %d", t);

        //---- compute Slater matrix elements in files
        for (int m=0; m<N; ++m) {
          for (int n=0; n<N; ++n) {
            corr[m][n][t] += sources[m]->Project(props[n]);
            corr2[m][n][t] += two_body.Run(props[m],props[n]);
          }
        }

        //---- Compute propagators on the current background configuration
        for (int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        for (int n=0; n<N; ++n) {
          //VRB.Debug(fname, "Propagator %d of %d\n", n, N);
          prop[n]->Run(hamiltonians);
        }

        //---- Generate a new field configuration
        lattice2.Refresh();

      }

    }

    //---- Store accumulated results in file 
    for (int m=0; m<N; ++m) {
      for (int n=0; n<N; ++n) {
        for(int t=0; t<GJP.Tsites(); ++t) {
          //---- Write single particle propagator to disk
          file[m][n] << setprecision(PREC) << corr[m][n][t].real()/(Float)unload_period << " ";
          file[m][n] << setprecision(PREC) << corr[m][n][t].imag()/(Float)unload_period << " ";
          corr[m][n][t]=0.0;
          //---- Write multi-particle propagator to disk
          file2[m][n] << setprecision(PREC) << corr2[m][n][t].real()/(Float)unload_period << " ";
          file2[m][n] << setprecision(PREC) << corr2[m][n][t].imag()/(Float)unload_period << " ";
          corr2[m][n][t]=0.0;
        }
        file[m][n] << endl;
        file2[m][n] << endl;
      }
    }

  }

  //---- Delete propagator memory allocation
  for (int n=0; n<N; ++n) {
    delete prop[n];
  }

  //---- Delete one-body memory allocation
  for (int n=0; n<N; ++n) {
    delete sources[n];
  }

  //---- Close data files 
  for (int m=0; m<N; ++m) {
    for (int n=0; n<N; ++n) {
      file[m][n].close();
      file2[m][n].close();
    }
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

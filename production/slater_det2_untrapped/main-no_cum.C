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
#include "slater_det2.h"
#include "verbose.h"
#include "error.h"
#include "momenta.h"
#include "global_job_parameter.h"
#include "dispersion.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "slater_det2_untrapped";

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

//  PotentialArg potential_arg;
//  potential_arg.Decode("args/potential.arg");

  KineticArg kinetic_arg;
  kinetic_arg.Decode("args/kinetic.arg");

  PropagatorArg prop_arg;
  prop_arg.Decode("args/propagator.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  MomentaArg momenta_arg;
  momenta_arg.Decode("args/momenta.arg");

  //---- Perform some checks

  if (one_body_arg.source_type != SOURCE_TYPE_MOM) {
    VRB.Warn(fname, "one_body_arg.source_type not supported:");
    VRB.Warn(fname, "\tusing SOURCE_TYPE_MOM instead.");
    one_body_arg.source_type = SOURCE_TYPE_MOM;
  }

  if (momenta_arg.dispersion_type != kinetic_arg.dispersion_type) {
    VRB.Warn(fname, "momenta_arg.dispersion_type and prop_arg.dispersion_type disagree:");
    VRB.Warn(fname, "\tSetting momenta_arg.dispersion_type to prop_arg.dispersion_type.");
    momenta_arg.dispersion_type = kinetic_arg.dispersion_type;
  }

  if (momenta_arg.mass != kinetic_arg.mass) {
    VRB.Warn(fname, "momenta_arg.mass and prop_arg.mass disagree:");
    VRB.Warn(fname, "\tSetting momenta_arg.mass to prop_arg.mass.");
    momenta_arg.mass = kinetic_arg.mass;
  }

  if (evo_arg.unload_period<1) {
    ERR.General(fname,"evo_arg.unload_period must be greater than zero.");
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);
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

  //---- Determine lattice momenta below the fermi surface 
  Momenta momenta(momenta_arg);
  int N = momenta.MomentaCount();
  vector<int> momentum;

  if (!Comms::Rank()) {

    Dispersion dispersion(kinetic_arg.dispersion_type, kinetic_arg.mass, CUTOFF_TYPE_HARD);

    //---- Save momenta information to file
    ofstream file;
    char f_name[99];
    sprintf(f_name, "%s", "results/sources"); 
    file.open(f_name,ios_base::app);
    for (int n=0; n<N; ++n) {
      momentum = momenta.GetMomentum(n);
      file << n << " (" << momenta.OppositeParityIndex(n) << "): ";
      file << momentum[0] << " " << momentum[1] << " " << momentum[2] << " ";
      file << log(1+dispersion.Get(momentum[0],momentum[1],momentum[2])) << "\n";
    }
    file.close();

  }

  //---- Some measurement classes for computing Slater determinants of one particle wave-functions
  OneBody* one_bodies[N];
  for (int n=0; n<N; ++n) {
    momentum = momenta.GetMomentum(n);
    one_bodies[n] = new OneBody(one_body_arg);
    one_bodies[n]->Set(momentum[0],momentum[1],momentum[2]);
  }

  TwoBody two_body(two_body_arg);
  two_body.Deform("args/two_body_def.arg");

  //---- Open files for measurments (files labeled by time separation; each file will contain data for the entire ensemble)
  FILE* efile[N];
  FILE* ofile[momenta.NumShells()][N];
  char f_name[99];
  for (int n=0; n<N; ++n) {
    sprintf(f_name, "%s.%d", "results/slater2_even", n+1);
    efile[n] = Fopen(IO_TYPE_ROOT,f_name,"a");
    for (int s=0; s<momenta.NumShells(); ++s) {
      if (n>=momenta.AccShellCount(s-1)) {
        sprintf(f_name, "%s_%d.%d", "results/slater2_odd", s, n+1); 
        ofile[s][n] = Fopen(IO_TYPE_ROOT,f_name,"a");
      }
    }
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr_even[N][GJP.Tsites()];
  complex<Float> corr_odd[momenta.NumShells()][N][GJP.Tsites()];
  for (int n=0; n<N; ++n)
  for (int t=0; t<GJP.Tsites(); ++t) {
    corr_even[n][t]=0.0;
    for (int s=0; s<momenta.NumShells(); ++s) {
      if (n>=momenta.AccShellCount(s-1)) {
        corr_odd[s][n][t]=0.0;
      }
    }
  }

  //---- Instantiate propagator classes
  prop_arg.n_fermions = N;
  Propagator propagator(prop_arg);

  //--- Set pointers to propagators
  Float* uprops[N];
  Float* dprops[N];
  for (int n=0; n<N; ++n) {
    uprops[n] = propagator.Get(n);
    dprops[momenta.OppositeParityIndex(n)] = propagator.Get(n);
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

//        VRB.Debug(fname, "t = %d", t);

        //---- Compute Slater determinants of various sub-matrices
        SlaterDet2 slater_det(uprops, dprops, one_bodies, NULL, &two_body, N);

        for (int n=0; n<N; ++n) {
          //---- Case 1: n_up = n+1 , n_down = n+1
          corr_even[n][t] += slater_det.Run(n+1, -1, -1);
          //---- Case 2: n_up = n+1, n_down = n
          for (int s=0; s<momenta.NumShells(); ++s) {
            if (n>=momenta.AccShellCount(s-1)) {
              corr_odd[s][n][t] += slater_det.Run(n+1, -1, momenta.AccShellCount(s-1));
            } 
          } 
        } 

        //---- Compute propagators on the current background configuration
        for (int p=0; p<hamiltonians.size(); ++p) { hamiltonians[p]->Initialize(); }
        propagator.Run(hamiltonians);

        //---- Generate a new field configuration
        lattice.Refresh();

      }

    }

    //---- Perform global sums
    Comms::GlobalSum( &(corr_even[0][0]), N * GJP.Tsites() );
    Comms::GlobalSum( &(corr_odd[0][0][0]), momenta.NumShells() * N * GJP.Tsites() );

    //---- Store accumulated results in file 
    Float block_size = Comms::Size() * unload_period;
    for (int n=0; n<N; ++n) {
      for(int t=0; t<GJP.Tsites(); ++t) {
        //---- Write data for even fermions to file
        {
          //---- Store accumulated results in file 
          Fprintf(efile[n], "%.*e %.*e ", PREC, corr_even[n][t].real() / block_size, PREC, corr_even[n][t].imag() / block_size);
         
          //---- Reset accumulated corr
          corr_even[n][t]=0.0;
        }
        //---- Write data for odd fermions to file
        for (int s=0; s<momenta.NumShells(); ++s) {
          if (n>=momenta.AccShellCount(s-1)) {
            //---- Store accumulated results in file 
            Fprintf(ofile[s][n], "%.*e %.*e ", PREC, corr_odd[s][n][t].real() / block_size, PREC, corr_odd[s][n][t].imag() / block_size);
            //---- Reset accumulated corr
            corr_odd[s][n][t]=0.0;
          }
        }
      }
      Fprintf(efile[n], "\n");
      //--- Perform an occasional flush; this is needed because otherwise if the program quits unexpectedly, all data will be lost
      if (j % FLUSHRATE == 0 ) {
        Fflush(efile[n]);
      }
      for (int s=0; s<momenta.NumShells(); ++s) {
        if (n>=momenta.AccShellCount(s-1)) {
          Fprintf(ofile[s][n], "\n");
          if (j % FLUSHRATE == 0 ) {
            Fflush(ofile[s][n]);
          }
        }
      }
    }

  }

  //---- Delete wave function body memory allocation
  for (int n=0; n<N; ++n) {
    delete one_bodies[n];
  }

  //---- Delete data files 
  for (int n=0; n<N; ++n) {
    Fclose(efile[n]);
    for (int s=0; s<momenta.NumShells(); ++s) {
      if (n>=momenta.AccShellCount(s-1)) {
        Fclose(ofile[s][n]);
      }
    }
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

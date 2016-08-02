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
//  int N = 11;
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
//  int bn = 4;
//  TwoBody* two_bodies[bn];
//  for (int i=0; i<bn; ++i){
//    two_body_arg.lambda = 0.4 + 0.1*i;
//    two_bodies[bn] = new TwoBody(two_body_arg);
//    two_bodies[bn]->Deform("arg/two_body_def/arg");
//  }

  //---- Open files for measurments (files labeled by time separation; each file will contain data for the entire ensemble)
  FILE* efile[N];
  FILE* ofile[momenta.NumShells()][N];
  FILE* ocumfile[momenta.NumShells()][N];
  char f_name[99];
  for (int n=0; n<N; ++n) {
    sprintf(f_name, "%s.%d", "results/slater2_even", n+1);
    efile[n] = Fopen(IO_TYPE_ROOT,f_name,"a");
    for (int s=0; s<momenta.NumShells(); ++s) {
      if (n>=momenta.AccShellCount(s-1)) {
        sprintf(f_name, "%s_%d.%d", "results/slater2_odd", s, n+1); 
        ofile[s][n] = Fopen(IO_TYPE_ROOT,f_name,"a");
        sprintf(f_name, "%s_%d.%d", "results/slater2_odd_cum", s, n+1); 
        ocumfile[s][n] = Fopen(IO_TYPE_ROOT,f_name,"a");
      }
    }
  }

  //---- Some variables before outer loop
  int start = evo_arg.start;
  int unload_period = evo_arg.unload_period;
  int configurations = evo_arg.configurations;
  complex<Float> corr_even[N][GJP.Tsites()];
  complex<Float> corr_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> corr_odd_cum[momenta.NumShells()][N][GJP.Tsites()];
  int count[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z1_even[N][GJP.Tsites()];
  complex<Float> z2_even[N][GJP.Tsites()];
  complex<Float> z3_even[N][GJP.Tsites()];
  complex<Float> z4_even[N][GJP.Tsites()];
  complex<Float> z5_even[N][GJP.Tsites()];
  complex<Float> z6_even[N][GJP.Tsites()];
  complex<Float> z7_even[N][GJP.Tsites()];
  complex<Float> z8_even[N][GJP.Tsites()];
  complex<Float> z1_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z2_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z3_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z4_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z5_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z6_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z7_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> z8_odd[momenta.NumShells()][N][GJP.Tsites()];
  complex<Float> tmp;
  Float tmp_log;

  for (int n=0; n<N; ++n)
  for (int t=0; t<GJP.Tsites(); ++t) {
    corr_even[n][t]=0.0;
    z1_even[n][t]=0.0;
    z2_even[n][t]=0.0;
    z3_even[n][t]=0.0;
    z4_even[n][t]=0.0;
    z5_even[n][t]=0.0;
    z6_even[n][t]=0.0;
    z7_even[n][t]=0.0;
    z8_even[n][t]=0.0;
    for (int s=0; s<momenta.NumShells(); ++s) {
      if (n>=momenta.AccShellCount(s-1)) {
        corr_odd[s][n][t]=0.0;
        corr_odd_cum[s][n][t]=0.0;
        z1_odd[s][n][t]=0.0;
        z2_odd[s][n][t]=0.0;
        z3_odd[s][n][t]=0.0;
        z4_odd[s][n][t]=0.0;
        z5_odd[s][n][t]=0.0;
        z6_odd[s][n][t]=0.0;
        z7_odd[s][n][t]=0.0;
        z8_odd[s][n][t]=0.0;
	count[s][n][t]=0;
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

    VRB.Flow(fname, "Configuration %d of %d (with an unload_period of %d).",j+1,configurations, unload_period);

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
    VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);

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
          tmp = slater_det.Run(n+1, -1, -1);
          tmp_log = log (tmp.real());
          corr_even[n][t] += tmp;
          z1_even[n][t] += tmp_log;
          z2_even[n][t] += tmp_log*tmp_log;
          z3_even[n][t] += tmp_log*tmp_log*tmp_log;
          z4_even[n][t] += tmp_log*tmp_log*tmp_log*tmp_log;
          z5_even[n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
          z6_even[n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
          z7_even[n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
          z8_even[n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;

          //---- Case 2: n_up = n+1, n_down = n
          for (int s=0; s<momenta.NumShells(); ++s) {
            if (n>=momenta.AccShellCount(s-1)) {
              corr_odd[s][n][t] += slater_det.Run(n+1, -1, momenta.AccShellCount(s-1));
              tmp = slater_det.Run(n+1, -1, momenta.AccShellCount(s-1));
//	      if (tmp_odd.real() < 0) {
//		cout << "\n C =" << tmp_odd.real() <<" is negative at N =" << 2*n+1 << "and t =" << t << "\n";
//	      }
	      if (tmp.real() >0 ) {
                corr_odd_cum[s][n][t] += tmp;
          	tmp_log = log (tmp.real());
	        z1_odd[s][n][t] += tmp_log;
       		z2_odd[s][n][t] += tmp_log*tmp_log;
        	z3_odd[s][n][t] += tmp_log*tmp_log*tmp_log;
        	z4_odd[s][n][t] += tmp_log*tmp_log*tmp_log*tmp_log;
        	z5_odd[s][n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
        	z6_odd[s][n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
         	z7_odd[s][n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
        	z8_odd[s][n][t] += tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log*tmp_log;
  	      }
	      else {
		count[s][n][t] += 1;
//		cout << "\n C =" << tmp_odd.real() <<" is negative at N =" << 2*n+1 << "and t =" << t << "\n";
	      }
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

    //---- Store accumulated results in file 
    complex<Float> corr_glb_sum;
    complex<Float> corr_glb_odd_sum;
    complex<Float> z1_sum;
    complex<Float> z2_sum;
    complex<Float> z3_sum;
    complex<Float> z4_sum;
    complex<Float> z5_sum;
    complex<Float> z6_sum;
    complex<Float> z7_sum;
    complex<Float> z8_sum;

    for (int n=0; n<N; ++n) {
      for(int t=0; t<GJP.Tsites(); ++t) {
        //---- Write data for even fermions to file
        {
          //---- Need to perform a global average across all nodes
          corr_glb_sum = Comms::GlobalSum(corr_even[n][t]/(Float)unload_period);
          corr_glb_sum /= Comms::Size();
          z1_sum = Comms::GlobalSum(z1_even[n][t]/(Float)unload_period);
          z1_sum /= Comms::Size();
          z2_sum = Comms::GlobalSum(z2_even[n][t]/(Float)unload_period);
          z2_sum /= Comms::Size();
          z3_sum = Comms::GlobalSum(z3_even[n][t]/(Float)unload_period);
          z3_sum /= Comms::Size();
          z4_sum = Comms::GlobalSum(z4_even[n][t]/(Float)unload_period);
          z4_sum /= Comms::Size();
          z5_sum = Comms::GlobalSum(z5_even[n][t]/(Float)unload_period);
          z5_sum /= Comms::Size();
          z6_sum = Comms::GlobalSum(z6_even[n][t]/(Float)unload_period);
          z6_sum /= Comms::Size();
          z7_sum = Comms::GlobalSum(z7_even[n][t]/(Float)unload_period);
          z7_sum /= Comms::Size();
          z8_sum = Comms::GlobalSum(z8_even[n][t]/(Float)unload_period);
          z8_sum /= Comms::Size();

          //---- Then write the result to a file
          Fprintf(efile[n], "%.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e ", PREC, corr_glb_sum.real(), PREC, corr_glb_sum.imag(), PREC, z1_sum.real(), PREC, z2_sum.real(), PREC, z3_sum.real(), PREC, z4_sum.real(), PREC, z5_sum.real(), PREC, z6_sum.real(), PREC, z7_sum.real(), PREC, z8_sum.real());
         
          //---- Reset accumulated corr
          corr_even[n][t]=0.0;
          z1_even[n][t]=0.0;
          z2_even[n][t]=0.0;
          z3_even[n][t]=0.0;
          z4_even[n][t]=0.0;
          z5_even[n][t]=0.0;
          z6_even[n][t]=0.0;
          z7_even[n][t]=0.0;
          z8_even[n][t]=0.0;
        }
        //---- Write data for odd fermions to file
        for (int s=0; s<momenta.NumShells(); ++s) {
          if (n>=momenta.AccShellCount(s-1)) {
            //---- Need to perform a global average across all nodes
	    int unload_period_odd = unload_period - count[s][n][t];
            corr_glb_sum = Comms::GlobalSum(corr_odd[s][n][t]/(Float)unload_period);
            corr_glb_sum /= Comms::Size();
            corr_glb_odd_sum = Comms::GlobalSum(corr_odd_cum[s][n][t]/(Float)unload_period_odd);
            corr_glb_odd_sum /= Comms::Size();
            z1_sum = Comms::GlobalSum(z1_odd[s][n][t]/(Float)unload_period_odd);
            z1_sum /= Comms::Size();
            z2_sum = Comms::GlobalSum(z2_odd[s][n][t]/(Float)unload_period_odd);
            z2_sum /= Comms::Size();
            z3_sum = Comms::GlobalSum(z3_odd[s][n][t]/(Float)unload_period_odd);
            z3_sum /= Comms::Size();
            z4_sum = Comms::GlobalSum(z4_odd[s][n][t]/(Float)unload_period_odd);
            z4_sum /= Comms::Size();
            z5_sum = Comms::GlobalSum(z5_odd[s][n][t]/(Float)unload_period_odd);
            z5_sum /= Comms::Size();
            z6_sum = Comms::GlobalSum(z6_odd[s][n][t]/(Float)unload_period_odd);
            z6_sum /= Comms::Size();
            z7_sum = Comms::GlobalSum(z7_odd[s][n][t]/(Float)unload_period_odd);
            z7_sum /= Comms::Size();
            z8_sum = Comms::GlobalSum(z8_odd[s][n][t]/(Float)unload_period_odd);
            z8_sum /= Comms::Size();

            //---- Then write the result to a file
            Fprintf(ofile[s][n], "%.*e %.*e ", PREC, corr_glb_sum.real(), PREC, corr_glb_sum.imag());
            Fprintf(ocumfile[s][n], "%.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e %.*e ", PREC, corr_glb_odd_sum.real(), PREC, corr_glb_odd_sum.imag(), PREC, z1_sum.real(), PREC, z2_sum.real(), PREC, z3_sum.real(), PREC, z4_sum.real(), PREC, z5_sum.real(), PREC, z6_sum.real(), PREC, z7_sum.real(), PREC, z8_sum.real(), PREC, (Float)unload_period_odd);
            //---- Reset accumulated corr
            corr_odd[s][n][t]=0.0;
            corr_odd_cum[s][n][t]=0.0;
            z1_odd[s][n][t]=0.0;
            z2_odd[s][n][t]=0.0;
            z3_odd[s][n][t]=0.0;
            z4_odd[s][n][t]=0.0;
            z5_odd[s][n][t]=0.0;
            z6_odd[s][n][t]=0.0;
            z7_odd[s][n][t]=0.0;
            z8_odd[s][n][t]=0.0;
            count[s][n][t]=0;
          }
        }
      }
      Fprintf(efile[n], "\n");
      if (j % FLUSHRATE == 0 ) { Fflush(efile[n]); }
      for (int s=0; s<momenta.NumShells(); ++s) {
        if (n>=momenta.AccShellCount(s-1)) {
          Fprintf(ofile[s][n], "\n");
          if (j % FLUSHRATE == 0 ) { Fflush(ofile[s][n]); }
          Fprintf(ocumfile[s][n], "\n");
          if (j % FLUSHRATE == 0 ) { Fflush(ocumfile[s][n]); }
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
        Fclose(ocumfile[s][n]);
      }
    }
  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

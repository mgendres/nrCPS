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
#include "config.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

void Print(complex<Float>* mat, int N)
{

  for (int i=0; i<N; ++i) {
    for (int j=0; j<N; ++j) {
      cout << mat[j+N*i] << " ";
    }
    cout << endl;
  }

}


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

  PropagatorArg prop_arg;
  prop_arg.Decode("args/propagator.arg");

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  //---- Perform some checks
  if (one_body_arg.source_type!=SOURCE_TYPE_SHO) {
    VRB.Warn(fname,"Setting one_body_arg.source_type to SOURCE_TYPE_SHO.");
    one_body_arg.source_type=SOURCE_TYPE_SHO;
  }

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);

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

  //---- Instantiate propagator classes
  Propagator propagator(prop_arg);

  //--- Set pointers to propagators at time slice t
  Float* props[N];
  for (int n=0; n<N; ++n) {
    props[n] = propagator.Get(n);
  }

  //---- Set propagator sources and inititialize
  for (int n=0; n<N; ++n) {
    propagator.Set(one_bodies[n]->Get(), n);
  }

  //---- Set up timing variables
  complex<Float>* slater;
  clock_t clock_base;
  Float  clock_elapsed;

#ifdef USE_GPU
  Float** dev_props; // Eventualy move this to constructor?
  Cuda::Malloc( (void**)&dev_props, N*sizeof(Float*) );
  Cuda::MemCopy(dev_props, props, N*sizeof(Float*), cudaMemcpyHostToDevice);
#endif 

  for (int n=1; n<=N; ++n) {

    slater = (complex<Float> *) malloc(n*n*sizeof(complex<Float>) );

    clock_base = clock();


    for(int t=0; t<GJP.Tsites(); ++t) {

#ifdef USE_GPU
      Cuda::InnerProduct( dev_props , dev_props , slater, n, GJP.Vol());
#else
      for (int i=0; i<n; ++i) {
        for (int j=0; j<n; ++j) {
          slater[j+n*i] = one_bodies[i]->Project(props[j] );
        }
      }
#endif 

    }

    clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;

    VRB.Flow(fname, "%d %f", n, clock_elapsed);
    //VRB.Flow(fname, "{%d, %f}, ", n, clock_elapsed);
    //Print(slater, n);

    free(slater);

  }

  //---- Delete one-body memory allocation
  for (int n=0; n<N; ++n) {
    delete one_bodies[n];
  }

#ifdef USE_GPU
  Cuda::Free(dev_props);
#endif 

  Comms::Finalize();

  return(EXIT_SUCCESS);
}



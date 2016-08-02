#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
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

int main()
{

  const char* fname = "int main()";
  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg;
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  PropagatorArg prop_arg;
  prop_arg.Decode("args/prop.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  if (prop_arg.n_fermions != 2) {ERR.General(fname, "Number of prop_arg.n_fermions must equal two."); }
  Propagator propagator(prop_arg);

  TwoBody two_body(two_body_arg);

  two_body.Deform("args/two_body_def.arg");

  clock_t clock_base = clock();
  Float  clock_elapsed;

  for (int i=0; i<GJP.Tsites(); i++) { two_body.Run( propagator.Get(0), propagator.Get(1) ); }

  clock_elapsed = Float(clock()-clock_base)/CLOCKS_PER_SEC;
  VRB.Flow(fname, "Time elapsed is %f seconds.", clock_elapsed);

  return(EXIT_SUCCESS);

}

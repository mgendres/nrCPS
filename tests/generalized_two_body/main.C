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

  GeneralizedTwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");

  GeneralizedTwoBody two_body(two_body_arg);

  return(EXIT_SUCCESS);

}

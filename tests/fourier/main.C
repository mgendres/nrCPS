#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "fourier.h"

int main()
{

  const char* fname = "int main()";
  cout << "Here we go!\n\n";

  //---- Import simulation parameters

  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");

  DoArg do_arg;
  do_arg.Decode("args/do.arg");

  //---- Perform initialization
  VRB.SetLevel(verbose_arg);
  GJP.Initialize(do_arg);

  int N = GJP.Tsites();

  Fourier fourier(N);
  fourier.Forward();
  fourier.Backward();

  return(EXIT_SUCCESS);

}

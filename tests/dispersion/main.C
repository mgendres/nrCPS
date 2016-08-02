#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
using namespace std;
#include "enum.h"
#include "arg.h"
#include "verbose.h"
#include "dispersion.h"
#include "global_job_parameter.h"


int main()
{

  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg; // Container for lattice parameters
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  MomentaArg momenta_arg; // container for Momenta parameters
  momenta_arg.Decode("args/momenta.arg");

  Dispersion dispersion(momenta_arg.dispersion_type, momenta_arg.mass, CUTOFF_TYPE_HARD);
  Dispersion dispersion2(momenta_arg.dispersion_type, momenta_arg.mass, CUTOFF_TYPE_NONE);

  for (int i=0; i<GJP.Vol(); ++i) {
    cout << dispersion.Get(i) << " ";
    cout << dispersion2.Get(i) << "\n";
  }


  return(EXIT_SUCCESS);
}

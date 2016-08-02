#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <vector>
using namespace std;
#include "constants.h"
#include "enum.h"
#include "arg.h"
#include "verbose.h"
#include "momenta.h"
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

  //---- Open files for measurments
  //---- row = config number, col = time separation
  ofstream file;
  char f_name[99];
  sprintf(f_name, "%s", "results/kappa_inv_sq"); 
  file.setf(ios_base::scientific);
  file.open(f_name,ios_base::app);


  Float xi;
  Float kappa_inv_sq;
  Float lambda;

  Float min = 0.01;
  Float max = 3.0;
  Float inc = 0.00005;
  for (lambda=min; lambda<max; lambda+=inc) {

    kappa_inv_sq = 0.0;
    for (int i=0; i<GJP.Vol(); ++i) {
      xi = 1.0 + dispersion.Get(i);
      kappa_inv_sq += 1.0/(lambda*xi*xi-1.0);
    }
    kappa_inv_sq /= GJP.Vol();

    file << setprecision(PREC) << -log(lambda) << " " << kappa_inv_sq << endl;
  }

  file.close();

  return(EXIT_SUCCESS);

}

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

  Momenta momenta(momenta_arg);
  Dispersion dispersion(momenta_arg.dispersion_type, momenta_arg.mass, CUTOFF_TYPE_HARD);
  vector<int> momentum;

//  cout << "Fermi energy: " << momenta_arg.fermi_energy << "\n";
//  cout << "Number of momenta below Fermi surface: " << momenta.MomentaCount() << "\n";

//  cout << "Total number of shells (up to Fermi surface): " << momenta.NumShells() << endl;
//  cout << "Number of momenta in each shell:" << endl;
//  for (int i=0; i<momenta.NumShells(); ++i) {
//    cout << momenta.ShellCount(i) << endl;
//  }

  for (int s=-1; s<momenta.NumShells(); ++s) {
    cout << momenta.AccShellCount(s) << endl;
  }



  for (int n=0; n<momenta.MomentaCount(); ++n) {
    momentum = momenta.GetMomentum(n);
    cout << momentum[0] << " ";
    cout << momentum[1] << " ";
    cout << momentum[2] << " ";
    cout << "    ";
    momentum = momenta.GetMomentum(momenta.OppositeParityIndex(n));
    cout << momentum[0] << " ";
    cout << momentum[1] << " ";
    cout << momentum[2] << " ";
    cout << "    ";

    cout << log( 1.0 + dispersion.Get(momentum[0],momentum[1],momentum[2])) << "\n";

  }

//  for (int n=0; n<GJP.Vol(); ++n) {
//    cout << log( 1.0 + dispersion.Get(n)) << " ";
//    cout << 1.0/(1 + dispersion.Get(n)) << "\n";
//  }

  return(EXIT_SUCCESS);
}

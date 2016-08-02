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
#include "random.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "lattice";

  VRB.Flow(fname, "Running production code %s, compiled on %s at %s. ", prog_name, __DATE__, __TIME__ );
  VRB.LibInfo(fname);

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg;
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  LatticeArg lattice_arg;
  lattice_arg.Decode("args/lattice.arg");

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  //---- Initialize random number generator
  Random rng(random_arg);

  //---- Instantiate the lattice (contains two and three-body fields, initialized when instantiated)
  Lattice lattice(lattice_arg, &rng);

  lattice.Refresh();
//  for (int i=0; i<lattice.FieldSize(); ++i) {
//    cout << lattice[i] << endl;
//    lattice[i]=0;
//  }

  cout << lattice.Mean() << endl;
  cout << lattice.Variance() << endl;

  lattice.Refresh(5.0);
//  for (int i=0; i<lattice.FieldSize(); ++i) {
//    cout << lattice[i] << endl;
//    lattice[i]=0;
//  }

  cout << lattice.Mean() << endl;
  cout << lattice.Variance() << endl;

  lattice.Refresh( complex<Float>(5.0, 5.0) );
//  for (int i=0; i<lattice.FieldSize(); ++i) {
//    cout << lattice[i] << endl;
//    lattice[i]=0;
//  }

  cout << lattice.Mean() << endl;
  cout << lattice.Variance() << endl;



//  for (int i=0; i<lattice.FieldSize(); ++i) {
//    cout << lattice[i] << endl;
//  }

//  cout << lattice.Mean() << endl;
//  cout << lattice.Variance() << endl;

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

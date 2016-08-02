#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
using namespace std;
#include "enum.h"
#include "arg.h"
#include "random.h"
#include "verbose.h"
#include "sysfunc.h"

int main()
{

  //---- Establish communications (needed to pass memory check)
  Comms::Initialize();

  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  cout << "Random Type: " <<  random_arg.random_type << "\n";
  cout << "Seed Type: " <<  random_arg.seed_type << "\n";
  cout << "Seed: " << random_arg.seed << "\n";

  //---- Initialize random number generator
  Random rng(random_arg);

  int N = 100000;

  if (0) {
    for (int i=0; i<N; ++i) { cout << rng.Uniform() << endl; }
  }

  if (0) {
    for (int i=0; i<N; ++i) { cout << rng.Normal(1.3, 0.7) << endl; }
  }

  if (1) {
    for (int i=0; i<N; ++i) { cout << rng.LogNormal(0.3, 0.7) << endl; }
  }

  //---- Close communications (needed to pass memory check)
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

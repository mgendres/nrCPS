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
 
  //---- Test output of a single stream
  if (0) {

      //---- Initialize random number generator
      Random rng(random_arg);
 
      cout << "Uniform (single stream):\n";
      for (int i=0; i<10; ++i) { cout << rng.Uniform() << endl; }
  
  }

  //---- Test independent streams 
  if (1) {
  
      //---- Initialize random number generator
      Random rng(random_arg);
      Random rng2(random_arg);
  
      cout << "\nUniform (two streams):\n";
      for (int i=0; i<10; ++i) { cout << rng.Uniform() << " " << rng2.Uniform() << endl; }
  
      cout << "\nGaussian (two streams):\n";
      for (int i=0; i<10; ++i) { cout << rng.Gauss(0.5) << " " << rng2.Gauss(0.5) << endl; }
  
      cout << "\nZ(4) (two streams):\n";
      for (int i=0; i<10; ++i) { cout << rng.Z(4) << " " << rng2.Z(4) << endl; }
  }

  //---- Test timing
  if (1) {

      //---- Initialize random number generator
      Random rng(random_arg);
      int n = 1000 ;
      cout << "Uniform (single stream, " << n << " calls):\n";
      for (int i=0; i<n; ++i) rng.Uniform();
  
  }

  //---- Close communications (needed to pass memory check)
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

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

  cout << "Here we go!\n\n";

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  RandomArg random_arg;
  random_arg.Decode("args/random.arg");

  {
    random_arg.seed_type = SEED_TYPE_INPUT;
    Random rng(random_arg);

    //---- Let the generator run a while
    for (int i=0; i<19; ++i) {
      rng.Uniform();
    }

    //----Write the state to file
    rng.WriteState(random_arg.file_stem.c_str());

    //---- Okay, here's the test:
    cout << endl;
    for (int i=0; i<11; ++i) {
      cout << rng.Uniform() << endl;
    }
    cout << endl;


  }

  {
    random_arg.seed_type = SEED_TYPE_FILE;
    Random rng(random_arg);

    //---- Okay, here's the test:
    cout << endl;
    for (int i=0; i<11; ++i) {
      cout << rng.Uniform() << endl;
    }
    cout << endl;


  }

  if (0) {
    Random rng(random_arg);
  
    //---- Let the generator run a while
    for (int i=0; i<19; ++i) {
      rng.Uniform();
    }
  
    //----Write the state to file
    rng.WriteState(random_arg.file_stem.c_str());
  
    //---- Okay, here's the test:
    cout << endl;
    for (int i=0; i<11; ++i) {
      cout << rng.Uniform() << endl;
    }
    cout << endl;
  
    //----Read the state to file
    rng.ReadState(random_arg.file_stem.c_str());
  
    //---- Okay, here's the test:
    cout << endl;
    for (int i=0; i<11; ++i) {
      cout << rng.Uniform() << endl;
    }
    cout << endl;

  }

  Comms::Finalize();

  return(EXIT_SUCCESS);
}

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
#include "sysio.h"

int main()
{

  //---- Establish communications (needed to pass memory check)
  Comms::Initialize();

  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  RandomArg random_arg("args/random.arg");

  cout << "Random Type: " <<  random_arg.random_type << "\n";
  cout << "Seed Type: " <<  random_arg.seed_type << "\n";
  cout << "Seed: " << random_arg.seed << "\n";
 
  //---- Initialize random number generator
  Random rng(random_arg);

  int T = 1000; 
  Float g = 0.5;
  int Nconf = 50000;
  int Nensemble = T;

  //---- Open files for measurments
  FILE* file;
  char f_name[99];
  sprintf(f_name, "%s", "results/data"); 
  file = Fopen(IO_TYPE_ROOT,f_name,"a");

  for (int i=0; i<Nensemble; ++i) {
 
    Float corr[T];
    for (int t=0; t<T; ++t) { corr[t]=0.0; }
    Float c;
 
    for (int n=0; n<Nconf; ++n) {
      c = 1.0;
      for (int t=0; t<T; ++t) {
        c *= (1.0 + g*(2.0*rng.Uniform() -1.0) );
        corr[t] += c;
      }
    }
 
    for (int t=0; t<T; ++t) { corr[t] /= Nconf; }
 
//    for (int t=0; t<T; ++t) { cout << corr[t] << endl; }
 
    for (int t=0; t<T; ++t) {
      Fprintf(file, "%.*e ", PREC, corr[t]);
      corr[t] = 0.0;
    }
    Fprintf(file, "\n");

  }

  //---- Close data file
  Fclose(file);

  //---- Close communications (needed to pass memory check)
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

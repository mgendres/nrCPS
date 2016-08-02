#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>

using namespace std;
#include "uniform_deviates.h"


void TestIO(UniformDeviate*);

int main()
{

  cout << "Here we go!\n\n";

  //long seed = time((time_t *)NULL);
  long seed = 32767;

  if (1) {
    cout << "ran0 with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ran0 rng;
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }

  if (1) {
    cout << "ran1 with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ran1 rng;
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }



  if (1) {
    cout << "ran2 with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ran2 rng;
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }

  if (1) {
    cout << "ran3 with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ran3 rng;
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }

  if (1) {
    cout << "ranlxs(0) with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ranlxs rng(0);
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }

  if (1) {
    cout << "ranlxd(1) with seed: " << seed << endl; 

    //---- Initialize random number generators
    Ranlxd rng(1);
    rng.Init(seed);
    cout << rng.TypeQ() << endl;

    //---- Run random number generators
    for (int i=0; i<10; ++i) {
      cout << rng.Run() << endl;
    }

    TestIO(&rng);

  }




  return(EXIT_SUCCESS);
}


void TestIO(UniformDeviate* rng)
{

  int n=313;
  Float x[n];
  Float y[n];
  int test=0;

  //---- Break in the generator
  for (int i=0; i<10000; ++i) {
    rng->Run();
  }

  //---- Test I/O
  long *state;
  state = (long *)malloc(rng->StateSize()*sizeof(long));
  rng->GetState(state);

  for (int k=0; k<100; ++k) {
    for (int j=0; j<n; ++j) {
      x[j] = rng->Run();
    }
  }

  rng->SetState(state);

  for (int k=0; k<100; ++k) {
    for (int j=0; j<n; ++j) {
      y[j] = rng->Run();
    }
  }
 
  for (int k=0; k<n; ++k)
  {
     if (x[k]!=y[k])
        test=1; 
  }

  if (test==1)
  {
    printf("\n");
    printf("I/O routines for random generator do not work properly\n");
    printf("=> do not use on this machine\n");
    printf("\n");
  } else {
    printf("\n");
    printf("I/O routines for random generator work properly\n");
    printf("\n");
  }

  free(state);

}

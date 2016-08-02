#include <stdlib.h>
#include <stdio.h>
#include <iostream>
#include <string>
#include <iomanip>
#include <math.h>
using namespace std;
#include "constants.h"
#include "arg.h"
#include "verbose.h"
#include "global_job_parameter.h"
#include "one_body.h"
#include "two_body.h"

int main()
{

  cout << "Here we go!\n\n";

  //---- Import simulation parameters
  VerboseArg verbose_arg;
  verbose_arg.Decode("args/verbose.arg");
  VRB.SetLevel(verbose_arg);

  DoArg do_arg;
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg); 

  OneBodyArg one_body_arg;
  one_body_arg.Decode("args/one_body.arg");

  TwoBodyArg two_body_arg;
  two_body_arg.Decode("args/two_body.arg");


  if (1) {

    cout << "MOM sources:" << endl;
    one_body_arg.source_type = SOURCE_TYPE_MOM;

    //int N = 7;
    int N = 6;
    OneBody* one_body[N];
    TwoBody two_body(two_body_arg);

    char f_name[99];
    for (int i=0; i<N; ++i) {
      one_body[i] = new OneBody(one_body_arg);
      sprintf(f_name, "%s.%d", "mom/source", i); 
      one_body[i]->Set(f_name);
    }
 
 
    cout << endl;
 
    cout << "One body::" << endl;
    for (int i=0; i<N; ++i) {
      for (int j=0; j<N; ++j) {
        cout << one_body[i]->Project(one_body[j]->Get()) << " ";
      }
      cout << endl;
    }
    cout << endl;

    cout << "Two body::" << endl;
    for (int i=0; i<N; ++i) {
      for (int j=0; j<N; ++j) {
        cout << two_body.Run(one_body[i]->Get(), one_body[j]->Get()) << " ";
      }
      cout << endl;
    }
    cout << endl;

    for (int i=0; i<N; ++i) {
      delete one_body[i];
    }

  }


  if ( !GJP.APBCQuery() ) {

    cout << "SHO sources:" << endl;
    one_body_arg.source_type = SOURCE_TYPE_SHO;

    int N = 8;
    OneBody* one_body[N];
    TwoBody two_body(two_body_arg);

    char f_name[99];
    for (int i=0; i<N; ++i) {
      one_body[i] = new OneBody(one_body_arg);
      sprintf(f_name, "%s.%d", "sho/source", i); 
      one_body[i]->Set(f_name);
    }
 
 
    cout << endl;
 
    for (int i=0; i<N; ++i) {
      for (int j=0; j<N; ++j) {
        cout << one_body[i]->Project(one_body[j]->Get()) << " ";
      }
      cout << endl;
    }
    cout << endl;

    cout << "Two body::" << endl;
    for (int i=0; i<N; ++i) {
      for (int j=0; j<N; ++j) {
        cout << two_body.Run(one_body[i]->Get(), one_body[j]->Get()) << " ";
      }
      cout << endl;
    }
    cout << endl;
 
    for (int i=0; i<N; ++i) {
      delete one_body[i];
    }

  } else {
    cout << "Cannot test SHO sources: APBCs are currently used." << endl;
  }


  return(EXIT_SUCCESS);

}

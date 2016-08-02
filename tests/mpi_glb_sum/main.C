#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "two_body";

  //---- Establish communications
  Comms::Initialize();

  int size = Comms::Size();
  int rank = Comms::Rank();

  //---- Print something
  cout << "\nI'm " << rank << " of " << size << "\n";


  complex<Float> z[5];
//  complex<Float> zgs[5];

  z[0] = complex<Float>(0.0,5.0);
  z[1] = complex<Float>(1.0,6.0);
  z[2] = complex<Float>(2.0,7.0);
  z[3] = complex<Float>(3.0,8.0);
  z[4] = complex<Float>(4.0,9.0);

//  Comms::GlobalSum(z,zgs,5);
  Comms::GlobalSum(z,5);

  if (Comms::Rank()==0) {
    for (int i=0; i<5; ++i) {
//      cout << zgs[i]/(Float)size << endl;
      cout << z[i]/(Float)size << endl;
    }
    cout << endl;
  }

  complex<Float> w[2][2][2];
//  complex<Float> wgs[2][2][2];

  w[0][0][0]= complex<Float>(0.0,8.0);
  w[0][0][1]= complex<Float>(1.0,9.0);
  w[0][1][0]= complex<Float>(2.0,10.0);
  w[0][1][1]= complex<Float>(3.0,11.0);
  w[1][0][0]= complex<Float>(4.0,12.0);
  w[1][0][1]= complex<Float>(5.0,13.0);
  w[1][1][0]= complex<Float>(6.0,14.0);
  w[1][1][1]= complex<Float>(7.0,15.0);

//  Comms::GlobalSum(&(w[0][0][0]), &(wgs[0][0][0]) ,8);
  Comms::GlobalSum(&(w[0][0][0]), 8);

  if (Comms::Rank()==0) {
    for (int i=0; i<8; ++i) {
//      cout << *( &wgs[0][0][0]+i) / (Float)size << endl;
      cout << *( &w[0][0][0]+i) / (Float)size << endl;

    }
    cout << endl;
  }

  //---- Close communications
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

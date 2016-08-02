#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "special_functions.h"
#include "constants.h"


int main()
{

  cout << "Here we go!\n\n";

  int n=10;
  Float x = 0.7;

  cout << "x = " << x << " and n = " << n << endl;

  cout << endl;
  cout << "method 1:" << endl;
  if (1) {
    Float h[n];
  
    HermiteH(h, n, x);
  
    for (int i=0; i<n; ++i) {
      cout << setprecision(PREC) << h[i] << endl;
    }
  }

  cout << endl;
  cout << "method 2:" << endl;
  if (1) {
    for (int i=0; i<n; ++i) {
      cout << setprecision(PREC) << HermiteH(i,x) << endl;
    }
  }


  return(EXIT_SUCCESS);
}

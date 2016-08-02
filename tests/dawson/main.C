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
  Float x = 0.1327;

  for (int i=-n-1; i<=n; ++i) {
    cout << setprecision(PREC) << i*x << " " << DawsonF(i*x) << endl;
  }

  return(EXIT_SUCCESS);
}

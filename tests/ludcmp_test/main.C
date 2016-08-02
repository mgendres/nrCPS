#include <stdlib.h>
#include <complex>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <math.h>
using namespace std;
#include "config.h"
#include "ludcmp.h"


int main()
{

  cout << "Here we go!\n\n";

  int n = 10;

  //---- Allocate memory for fermion matrix
  complex<Float> mat[n][n];

  cout << "Computing the determinent of the following matrix:" << endl;
  cout << endl;
  //---- Construct fermion matrix
  for (int s=0; s<n; ++s) {
    for (int t=0; t<n; ++t) {
      mat[s][t] = complex<Float>(1.0/(s+t+2.0),0.0);
      printf("%f ",real(mat[s][t]));
      //printf("%Lf ",real(mat[s][t]));
    }
    cout << endl;
  }
  cout << endl;

  //---- Compute LU decomposition
  int indx[n];
  Float d;
  ludcmp(&mat[0][0],n,indx,&d);


  //---- Compute Slater det from LU decomposition
  if (1) {
    complex<Float> det(d,0.0);
    for(int i=0; i<n; ++i) {
      det *= mat[i][i];
    }
    //printf("det = (%Le, %Le)\n",real(det),imag(det));
    printf("det = (%e, %e)\n",real(det),imag(det));
  }

  cout << "det = 1.17137155298420180 * 10^-58 (Mathematica, e^18 precision)\n";

  //---- Compute log of Slater det from LU decomposition
  if (1) {
    complex<Float> ldet;
    ldet = log (complex<Float> (d,0.0));
    for(int i=0; i<n; ++i) {
      ldet += log(mat[i][i]);
    }
    printf("log(det) = (%e, %e)\n",real(ldet),imag(ldet));
    //printf("log(det) = (%Le, %Le)\n",real(ldet),imag(ldet));
  }

  cout << "ldet = -133.391760063906166 (Mathematica, e^18 precision)\n";

  return(EXIT_SUCCESS);

}

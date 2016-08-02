#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "slater_det2.h"
#include "ludcmp.h"
#include "verbose.h"
#include "error.h"
#include "enum.h"

void SlaterDet2::Initialize(int n)
{

  N = n;

  //---- Initialize pointers
  psi1 = NULL;
  psi2 = NULL;

  //---- Allocate memory for fermion matrix
  phi = (complex<Float> *) malloc(N*N*sizeof(complex<Float>)); //Allocate memory for field

  //---- Default determinent setting
  SetDetType(DET_TYPE_STANDARD);

}

SlaterDet2::SlaterDet2(Float** props1, Float** props2, OneBody** one_body1, OneBody** one_body2, TwoBody* two_body, int n)
// Note that one_body1 and/or one_body2 may be NULL pointers.
{

  Initialize(n);

  if (one_body1 !=NULL) {
    psi1 = (complex<Float> *) malloc(N*N*sizeof(complex<Float>)); //Allocate memory for field
    for (int s=0; s<N; ++s) {
      for (int t=0; t<N; ++t) {
        psi1[t+N*s] = one_body1[s]->Project(props1[t]);
      }
    }
  }

  if (one_body2 !=NULL) {
    psi2 = (complex<Float> *) malloc(N*N*sizeof(complex<Float>)); //Allocate memory for field
    for (int s=0; s<N; ++s) {
      for (int t=0; t<N; ++t) {
        psi2[t+N*s] = one_body2[s]->Project(props2[t]);
      }
    }
  }

  for (int s=0; s<N; ++s) {
    for (int t=0; t<N; ++t) {
      phi[t+N*s] = two_body->Run(props1[s], props2[t]);
    }
  }

}


SlaterDet2::~SlaterDet2()
{

  free(phi);

  if (psi1!=NULL) {
    free(psi1);
  }

  if (psi2!=NULL) {
    free(psi2);
  }

}

void SlaterDet2::SetDetType(DetType type) { det_type = type; }

complex<Float> SlaterDet2::Run(int n, int s0, int t0)
// s0 and t0 specify a specific row or column to strike out.
// Only s0 or t0 may be non-negative, but not both.
// If s0 is non-negative, then the type 1 fermion with quantum number s0 is removed.
// Similarly, if t0 is non-negative, then the type 2fermion with quantum number t0 is removed.
// This function will eventually be generalized to a list of s0 and t0s.
{

  char* fname = "complex<Float> SlaterDet2::Run(int, int, int)";
  VRB.Func(fname);

  //---- Perform lots of checks
  if (n>N) { ERR.General(fname, "Sub-matrix exceeds size of matrix."); }
  if (n<1) { ERR.General(fname, "Specified size of sub-matrix is not positive."); }
  if ((s0>=0)&&(t0>=0)) { ERR.General(fname, "s0 and t0 cannot both be non-negative."); }

  //---- Construct submatrix
  complex<Float>** mat = new complex<Float>* [n];
  for(int i = 0; i < n ; ++i)
    mat[i] = new complex<Float> [n];

  if ((t0>=0)&&(t0<n)) {
    if (psi1==NULL) { ERR.General(fname, "t0 cannot be non-negative; psi1==NULL."); }
    for (int s=0; s<n; ++s) {
      mat[s][t0] = psi1[s+t0*N];
    }
  }

  if ((s0>=0)&&(s0<n)) {
    if (psi2==NULL) { ERR.General(fname, "s0 cannot be non-negative; psi2==NULL."); }
    for (int t=0; t<n; ++t) {
      mat[s0][t] = psi2[t+s0*N];
    }
  }

  for (int s=0; s<n; ++s) {
    for (int t=0; t<n; ++t) {
      if ( (s!=s0)&&(t!=t0)) { mat[s][t] = phi[t+N*s]; }
    }
  }

  //---- Compute LU decomposition
  int indx[n];
  Float d;
  ludcmp(&mat[0][0],n,indx,&d);

  //---- Compute Slater det from LU decomposition
  if (det_type == DET_TYPE_STANDARD) {
    complex<Float> det(d,0.0);
    for(int i=0; i<n; ++i) {
      det *= mat[i][i];
    }
    return det;
  }


  //---- Compute log of Slater det from LU decomposition
  if (det_type == DET_TYPE_LOG) {
    complex<Float> ldet;
    ldet = log (complex<Float> (d,0.0));
    for(int i=0; i<n; ++i) {
      ldet += log(mat[i][i]);
    }
    return ldet;
  }

  for(int i = 0; i < n ; ++i)
    delete [] mat[i];
  delete [] mat;

}

complex<Float> SlaterDet2::SlaterMat(int s, int t) { return phi[t+N*s]; }


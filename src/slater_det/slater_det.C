#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;
#include "slater_det.h"
#include "ludcmp.h"
#include "verbose.h"
#include "error.h"
#include "enum.h"

void SlaterDet::Initialize(int n)
{

  mat_size = n;

  //---- Allocate memory for fermion matrix
  slater_mat = (complex<Float> *) malloc(mat_size*mat_size*sizeof(complex<Float>)); //Allocate memory for field

  //---- Default determinent setting
  SetDetType(DET_TYPE_STANDARD);

}

SlaterDet::SlaterDet(Float** props, OneBody** one_body, int n)
{

  Initialize(n);

  //---- Construct the Slater matrix M(i,j) = <i| e^{-TH} |j> out of single particle states
  //----  <i| pointed to by one_body[i] and  e^{-TH} |j> pointed to by props[i].

  for (int s=0; s<mat_size; ++s) {
    for (int t=0; t<mat_size; ++t) {
      slater_mat[t+mat_size*s] = one_body[s]->Project(props[t]);
    }
  }

}

SlaterDet::SlaterDet(Float** props1, Float** props2, TwoBody* two_body, int n)
{

  Initialize(n);

  //---- Construct the Slater matrix M(i,j) = <psi| e^{-TH} |ij> using the two particle
  //---- sink <psi| pointed to by two_body. The state e^{-TH}|ij> is decomposed into single
  //---- particle states e^{-TH}|i> pointed to by props1[i] and e^{-TH}|j> pointed to by
  //---- props2[j].

  for (int s=0; s<mat_size; ++s) {
    for (int t=0; t<mat_size; ++t) {
      slater_mat[t+mat_size*s] = two_body->Run(props1[s], props2[t]);
    }
  }

}

complex<Float> SlaterDet::Run(int n)
{

  const char* fname = "complex<Float> SlaterDet::Run(int)";
  VRB.Func(fname);

  if (n > mat_size) {
    ERR.General(fname, "Submatrix exceeds size of matrix.");
  }

  if (n<1) {
    VRB.Warn(fname, "n < 1, returning unity.");
    return 1.0;
  }

  //---- Construct Slater sub-matrix
  //---- The sub-matrix is given by the nxn corner of the full mat_size x mat_size matrix
//  complex<Float> mat[n][n];
  complex<Float>** mat = new complex<Float>* [n];
  for(int i = 0; i < n ; ++i)
    mat[i] = new complex<Float> [n];

  for (int s=0; s<n; ++s) {
    for (int t=0; t<n; ++t) {
      mat[s][t] = SlaterMat(s,t);
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


SlaterDet::~SlaterDet()
{
  free(slater_mat);
}

void SlaterDet::SetDetType(DetType type)
{
  det_type = type;
}

complex<Float> SlaterDet::SlaterMat(int s, int t) { return slater_mat[t+mat_size*s]; }


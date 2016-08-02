#include <iostream>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "slater_det3.h"
#include "ludcmp.h"
#include "verbose.h"
#include "error.h"
#include "enum.h"
#include "config.h"
#include "global_job_parameter.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#define NCUT 5
#endif

SlaterDet3::SlaterDet3(Float** _props1, Float** _props2, OneBody** _one_body1, OneBody** _one_body2, int n)
{

  N = n;

  //---- Default determinent setting
  SetDetType(DET_TYPE_STANDARD);

  //---- Initialize pointers
  phi = NULL;
  psi1 = NULL;
  psi2 = NULL;
  props1 = NULL;
  props2 = NULL;
  one_body1 = NULL;
  one_body2 = NULL;

  //---- Allocate memory for fermion matrices
  phi = (complex<Float> *) malloc(N*N*sizeof(complex<Float>));
  if (_one_body1 !=NULL) { psi1 = (complex<Float> *) malloc(N*N*sizeof(complex<Float>)); }
  if (_one_body2 !=NULL) { psi2 = (complex<Float> *) malloc(N*N*sizeof(complex<Float>)); }

  //---- Allocate pointers to fermion fields (propagators and sources)
#ifdef USE_GPU

  Cuda::Malloc( (void**)&props1, n*sizeof(Float*) );
  Cuda::MemCopy(props1, _props1, n*sizeof(Float*), cudaMemcpyHostToDevice);
  if (_one_body1!=NULL) {
    Cuda::Malloc( (void**)&one_body1, n*sizeof(Float*) );
    Float* p_address;
    for (int i=0; i<n; ++i) {
      p_address = _one_body1[i]->Get();
      Cuda::MemCopy( one_body1+i, &p_address, sizeof(Float*), cudaMemcpyHostToDevice);
    }
  }
 
  Cuda::Malloc( (void**)&props2, n*sizeof(Float*) );
  Cuda::MemCopy(props2, _props2, n*sizeof(Float*), cudaMemcpyHostToDevice);
  if (_one_body2!=NULL) {
    Cuda::Malloc( (void**)&one_body2, n*sizeof(Float*) );
    Float* p_address;
    for (int i=0; i<n; ++i) {
      p_address = _one_body2[i]->Get();
      Cuda::MemCopy( one_body2+i, &p_address, sizeof(Float*), cudaMemcpyHostToDevice);
    }
  }

#else

  props1 = _props1;
  if (_one_body1!=NULL) { one_body1 = _one_body1; }

  props2 = _props2;
  if (_one_body2!=NULL) { one_body2 = _one_body2; }

#endif

}


SlaterDet3::~SlaterDet3()
{

  if (phi!=NULL) { free(phi); }
  if (psi1!=NULL) { free(psi1); }
  if (psi2!=NULL) { free(psi2); }

#ifdef USE_GPU
  if (props1!=NULL) { Cuda::Free(props1); }
  if (props2!=NULL) { Cuda::Free(props2); }
  if (one_body1!=NULL) { Cuda::Free(one_body1); }
  if (one_body2!=NULL) { Cuda::Free(one_body2); }
#endif

}


void SlaterDet3::Initialize(TwoBody* two_body )
// Note that one_body1 and/or one_body2 may be NULL pointers.
{


  if (one_body1 !=NULL) {
#ifdef USE_GPU
    Cuda::InnerProduct(one_body1, props1, psi1, N, GJP.Vol() );
#else
    for (int s=0; s<N; ++s) {
      for (int t=0; t<N; ++t) {
        psi1[t+N*s] = one_body1[s]->Project(props1[t]);
      }
    }
#endif
  }

  if (one_body2 !=NULL) {
#ifdef USE_GPU
    Cuda::InnerProduct(one_body2, props2, psi2, N, GJP.Vol() );
#else
    for (int s=0; s<N; ++s) {
      for (int t=0; t<N; ++t) {
        psi2[t+N*s] = one_body2[s]->Project(props2[t]);
      }
    }
#endif
  }


#ifdef USE_GPU
  two_body->Run(props1, props2, phi, N); 
#else
  for (int s=0; s<N; ++s) {
    for (int t=0; t<N; ++t) {
      phi[t+N*s] = two_body->Run(props1[s], props2[t]);
    }
  }
#endif

}


void SlaterDet3::SetDetType(DetType type) { det_type = type; }

complex<Float> SlaterDet3::Run(int n, int s0, int t0)
// s0 and t0 specify a specific row or column to strike out.
// Only s0 or t0 may be non-negative, but not both.
// If s0 is non-negative, then the type 1 fermion with quantum number s0 is removed.
// Similarly, if t0 is non-negative, then the type 2fermion with quantum number t0 is removed.
// This function will eventually be generalized to a list of s0 and t0s.
{

  const char* fname = "complex<Float> SlaterDet3::Run(int, int, int)";
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

complex<Float> SlaterDet3::SlaterMat(int s, int t) { return phi[t+N*s]; }


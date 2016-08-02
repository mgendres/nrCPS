#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
using namespace std;
#include "arg.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"
#include "dispersion.h"
#include "special_functions.h"


complex<Float> TwoBody::Run(Float* prop1, Float* prop2)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*)";
//  VRB.Func(fname);

  // Compute sum_p psi(p) G(p) G(-p)
  int j;
  int k;
  Float re=0.0;
  Float im=0.0;


  for (int i=0; i<vol; ++i) {

//    cout << i << " " << wavefunc[i] << endl;

    j = 2*i;
    k = 2*opposite_parity_index[i];

    re += (prop1[j]*prop2[k] - prop1[j+1]*prop2[k+1]) * wavefunc[i];
    im += (prop1[j]*prop2[k+1] + prop1[j+1]*prop2[k]) * wavefunc[i];

  }

  return complex<Float>(re,im);

}

complex<Float> GeneralizedTwoBody::Run(Float* prop1, Float* prop2, OctahedralRep R)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*, OctahedralRep, const Vector3&)";
//  VRB.Func(fname);

  // Compute sum_q \pi_R(q) sum_p psi(p) G(p+q) G(-p)
  // Where the sum on q is over q_k = G(k).distinct_vector and G(k) are the 48 group elements of Oh


  Float re=0.0;
  Float im=0.0;

  int j, k;
  Float weight;

  for (int l=0; l<oh_group.GetSubshellDimension(); ++l) {
    weight = oh_group.GetWeight(R, l); 
    for (int i=0; i<GJP.Vol(); ++i) {
      j = 2 * shift_indices[l*vol+i];
      k = 2 * opposite_parity_index[i];

      re += weight * (prop1[j]*prop2[k] - prop1[j+1]*prop2[k+1]) * half_site_wavefunc[l*vol+i];
      im += weight * (prop1[j]*prop2[k+1] + prop1[j+1]*prop2[k]) * half_site_wavefunc[l*vol+i];

    } 
  }

  return complex<Float>(re,im);

}

complex<Float> GeneralizedTwoBody::Run2(Float* prop1, Float* prop2, OctahedralRep R)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*, OctahedralRep, const Vector3&)";
//  VRB.Func(fname);

  // Compute sum_q \pi_R(q) sum_p G(p + Q) G(-p)
  // Where the sum on q is over q_k = G(k).distinct_vector and G(k) are the 48 group elements of Oh
  // Ideally choose q to be a multiple of 2...

  Float re=0.0;
  Float im=0.0;

  int j, k;
  Float weight;
  Vector3 vec;

  for (int l=0; l<oh_group.GetSubshellDimension(); ++l) {
    weight = oh_group.GetWeight(R, l); 
    vec = oh_group.GetVector(l);
    j = 2 * Vector3ToIndex( vec );
    k = 2 * Vector3ToIndex( -1*vec );
    re += weight * (prop1[j]*prop2[k] - prop1[j+1]*prop2[k+1]);
    im += weight * (prop1[j]*prop2[k+1] + prop1[j+1]*prop2[k]);
  }

  return complex<Float>(re,im);

}


complex<Float> GeneralizedTwoBody2::Run(Float* prop1, Float* prop2)
{

//  const char* fname = "complex<Float> TwoBody::Run(Float*, Float*, OctahedralRep, const Vector3&)";
//  VRB.Func(fname);


  Float re=0.0;
  Float im=0.0;

  int j, k;

  for (int i=0; i<GJP.Vol(); ++i) {
    j = 2 * shift_index[i];
    k = 2 * opposite_parity_index[i];

    re += (prop1[j]*prop2[k] - prop1[j+1]*prop2[k+1]) * wavefunc[i];
    im += (prop1[j]*prop2[k+1] + prop1[j+1]*prop2[k]) * wavefunc[i];

  } 

  return complex<Float>(re,im);

}



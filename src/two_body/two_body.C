#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "arg.h"
#include "two_body.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"
#include "dispersion.h"
#include "special_functions.h"
#include "vector3.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

TwoBody::TwoBody(const TwoBodyArg &two_body_arg)
: vol( GJP.Vol() ),
  dims( GJP.GetDimensions() ),
  bcs( GJP.GetBoundaryConditions() )
{
  const char* fname = "TwoBody::TwoBody(const TwoBodyArg &)";
  VRB.Func(fname);

  wavefunc = (Float *) malloc(vol*sizeof(Float));
  opposite_parity_index= (int *) malloc(vol*sizeof(int));

  for (int i=0; i<vol; ++i) {
    opposite_parity_index[i] = Vector3ToIndex(-1* (IndexToVector3(i)+bcs) , dims);
  }

  //---- Construct two body wave function
  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_NONE) {
    for (int i=0; i<vol; ++i) { wavefunc[i] = 0.0; }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_UNIFORM) {
    for (int i=0; i<vol; ++i) { wavefunc[i] = 1.0/vol; }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_GND) {

    Dispersion dispersion1(two_body_arg.dispersion_type1, two_body_arg.mass1, CUTOFF_TYPE_HARD);
    Dispersion dispersion2(two_body_arg.dispersion_type2, two_body_arg.mass2, CUTOFF_TYPE_HARD);

    Float xi1;
    Float xi2;
    Float lambda = two_body_arg.lambda;

    for (int i=0; i<vol; ++i) {
      xi1 = 1.0 + dispersion1.Get(i);
      xi2 = 1.0 + dispersion2.Get(i);
 
      // Here, I use the fermion wave function, as derived in my note.
      // Note that psi(p) = xi(p)/(lambda*xi(p)^2-1) is the eigenstate of the transfer matrix.
      // The correlation function, however, is given by <final| D^{-1/2} T^N D^{-1/2} |initial>
      // Hence <p|final> = xi(p) psi(p), where xi is the same as D^{1/2} in momentum space).
      //
      // wavefunc[i] = xi^2/(lambda*xi^2-1.0);
      wavefunc[i] = 1.0;                     // This way avoids NAN when dispersion type is QUADRATIC, or PERFECT
      wavefunc[i] /= lambda - 1.0/(xi1*xi2); // so do it this way instead:
    }
  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_TGND) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    Float b = two_body_arg.lambda;
    Float psq;
    Float xi;

    for (int i=0; i<vol; ++i) {
      psq = 2.0*dispersion.Get(i); 

      if ( psq > GJP.CutoffSq() ) {

        wavefunc[i] = 0.0;

      } else {

        xi = sqrt(psq)/(2.0*b);

        if (xi < 1e-15) { // Wave function is non-singular at p=0, but need to treat this case separately
          wavefunc[i] = 1.0;
        } else {
          wavefunc[i] = DawsonF(xi)/xi;
        }

      }

    }

  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_PAIR1) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    Float bsq = two_body_arg.lambda;
    Float psq;

    for (int i=0; i<vol; ++i) {
      psq = 2.0*dispersion.Get(i); 
      wavefunc[i] = 1.0/(psq+bsq);
    }

  }

  if (two_body_arg.wavefunc_type==WAVEFUNC_TYPE_PAIR2) {

    Dispersion dispersion(DISPERSION_TYPE_QUADRATIC, 1.0, CUTOFF_TYPE_HARD);
    Float b = two_body_arg.lambda;
    Float psq;

    for (int i=0; i<vol; ++i) {
      psq = 2.0*dispersion.Get(i); 
      if (psq < 1e-15) {
        wavefunc[i] = 0.0; // Omit divergent contribution to wave funcion--must be treated separately
      } else {
        wavefunc[i] = exp(- b*sqrt(psq) )/psq;
      }
    }

  }

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_wavefunc, vol*sizeof(Float) );
  Cuda::Malloc( (void**)&dev_opposite_parity_index, vol*sizeof(int) );
  Cuda::MemCopy(dev_wavefunc, wavefunc, vol*sizeof(Float), cudaMemcpyHostToDevice);
  Cuda::MemCopy(dev_opposite_parity_index, opposite_parity_index, vol*sizeof(int), cudaMemcpyHostToDevice);

  //psum = (Float *) malloc(2*BLOCKS*sizeof(Float));
  Cuda::HostAlloc( (void**)&psum, 2*BLOCKS*sizeof(Float), cudaHostAllocDefault );
  Cuda::Malloc( (void**)&dev_psum, 2*BLOCKS*sizeof(Float) );

#endif

}

TwoBody::~TwoBody()
{

  const char* fname = "TwoBody::~TwoBody()";
  VRB.Func(fname);

#ifdef USE_GPU
  //free(psum);
  Cuda::FreeHost(psum);
  Cuda::Free(dev_psum);

  Cuda::Free( dev_wavefunc );
  Cuda::Free( dev_opposite_parity_index );
#endif

  free(wavefunc);
  free(opposite_parity_index);
}

void TwoBody::Deform(string vml_filename)
{

  const char* fname = "void TwoBody::Deform(string)";
  VRB.Func(fname);

  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) {
    ERR.FileR(fname,vml_filename.c_str());
  }

  int coord1;
  int coord2;
  int coord3;

  Float re;

  int j;

  //---- Need to add error checking, etc....
  while (true) {

    file >> coord1;
    file >> coord2;
    file >> coord3;

    j = Vector3ToIndex( Vector3(coord1, coord2, coord3) , dims);

    file >> re; 

    VRB.Flow(fname,"Appending %f to %d component of wave function.", re, j);
    wavefunc[j] += re;

    if ( file.eof() ) break; // Should be moved two lines up--if EOF, then accessing random memory.

  }

  file.close();

#ifdef USE_GPU
  Cuda::MemCopy(dev_wavefunc, wavefunc, vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif


}

GeneralizedTwoBody::GeneralizedTwoBody(const GeneralizedTwoBodyArg &two_body_arg)
: oh_group(two_body_arg.distinct_vec, 1),
  vol( GJP.Vol() ),
  dims( GJP.GetDimensions() ),
  bcs( GJP.GetBoundaryConditions() )
{
  const char* fname = "GeneralizedTwoBody::GeneralizedTwoBody(const TwoBodyArg &, const Vector3&)";
  VRB.Func(fname);

  int subshell_dimension = oh_group.GetSubshellDimension();
  opposite_parity_index= (int *) malloc(vol*sizeof(int));
  shift_indices = (int *) malloc(subshell_dimension*vol*sizeof(int));
  half_site_wavefunc = (Float *) malloc(subshell_dimension*vol*sizeof(Float));

  Vector3 site, shift, Q;
  Float psq;

  for (int i=0; i<vol; ++i) {
    site = IndexToVector3(i, dims);
    opposite_parity_index[i] = Vector3ToIndex(-1* (site+bcs) , dims);
    for (int l=0; l<subshell_dimension; ++l) {
      shift = oh_group.GetVector(l);
      shift_indices[l*vol+i] = Vector3ToIndex( site + shift , dims);

      // This mess is for evaluating the wavefunction
      Q = 2*site + shift + bcs;; // Wave function should be evaluated at Q/2 = site + shift/2
      Q = IndexToVector3( Vector3ToIndex( Q, 2*dims), 2*dims); // Move vector into BZ of doubled lattice
      Q = dims - abs(Q - dims);
      Q *= vol/dims;

      psq = Q.L2NormSquared();
      psq *= (PI/vol)*(PI/vol);


      //cout << "(" << site.GetX() << " " << site.GetY() << " " << site.GetZ() << ") : ";
      //cout << "(" << bcs.GetX() << " " << bcs.GetY() << " " << bcs.GetZ() << ") : ";
      //cout << "(" << shift.GetX() << " " << shift.GetY() << " " << shift.GetZ() << ") : ";
      //cout << psq << endl;

      switch (two_body_arg.wavefunc_type) {
        case WAVEFUNC_TYPE_NONE:
          half_site_wavefunc[l*vol+i] = 0.0;
          break;
        case WAVEFUNC_TYPE_UNIFORM:
          half_site_wavefunc[l*vol+i] = 1.0/vol;
          break;
        case WAVEFUNC_TYPE_GND:
          ERR.NotImplemented(fname, "WavefuncType not supported.");
          break;
        case WAVEFUNC_TYPE_TGND:
          if (psq < 1e-5) {
            half_site_wavefunc[l*vol+i] = 1.0;
          } else {
            psq = sqrt(psq)/(2.0* two_body_arg.lambda);
            half_site_wavefunc[l*vol+i] = DawsonF(psq)/psq;
          }
          break;
        case WAVEFUNC_TYPE_PAIR1:
          half_site_wavefunc[l*vol+i] = 1.0/(psq+two_body_arg.lambda);
          break;
        case WAVEFUNC_TYPE_PAIR2:
          if (psq < 1e-15) {
            half_site_wavefunc[l*vol+i] = 100.0;
          } else {
            half_site_wavefunc[l*vol+i] = exp(- two_body_arg.lambda*sqrt(psq) )/psq;
          }
          break;
        default:
          ERR.NotImplemented(fname, "WavefuncType not supported.");
      }

      // End mess

      //cout << half_site_wavefunc[l*vol+i] << endl;

    }

  }

}

GeneralizedTwoBody::~GeneralizedTwoBody()
{

  const char* fname = "GenerlizedTwoBody::~GenerlizedTwoBody()";
  VRB.Func(fname);

  free(half_site_wavefunc);
  free(opposite_parity_index);
  free(shift_indices);
}







GeneralizedTwoBody2::GeneralizedTwoBody2(const GeneralizedTwoBodyArg &two_body_arg)
: vol( GJP.Vol() ),
  dims( GJP.GetDimensions() ),
  bcs( GJP.GetBoundaryConditions() ),
  shift(two_body_arg.distinct_vec)
{
  const char* fname = "GeneralizedTwoBody2::GeneralizedTwoBody2(const TwoBodyArg &, const Vector3&)";
  VRB.Func(fname);

  opposite_parity_index= (int *) malloc(vol*sizeof(int));
  shift_index = (int *) malloc(vol*sizeof(int));
  wavefunc = (Float *) malloc(vol*sizeof(Float));

  Vector3 site, Q;
  Float psq;

  for (int i=0; i<vol; ++i) {
    site = IndexToVector3(i, dims);
    opposite_parity_index[i] = Vector3ToIndex(-1* (site+bcs) , dims);
    shift_index[i] = Vector3ToIndex( site + shift , dims);

    // This mess is for evaluating the wavefunction
    Q = 2*site + shift + bcs;; // Wave function should be evaluated at Q/2 = site + shift/2
    Q = IndexToVector3( Vector3ToIndex( Q, 2*dims), 2*dims); // Move vector into BZ of doubled lattice
    Q = dims - abs(Q - dims);
    Q *= vol/dims;

    psq = Q.L2NormSquared();
    psq *= (PI/vol)*(PI/vol);

    switch (two_body_arg.wavefunc_type) {
      case WAVEFUNC_TYPE_NONE:
        wavefunc[i] = 0.0;
        break;
      case WAVEFUNC_TYPE_UNIFORM:
        wavefunc[i] = 1.0/vol;
        break;
      case WAVEFUNC_TYPE_GND:
        ERR.NotImplemented(fname, "WavefuncType not supported.");
        break;
      case WAVEFUNC_TYPE_TGND:
        if (psq < 1e-5) {
          wavefunc[i] = 1.0;
        } else {
          psq = sqrt(psq)/(2.0* two_body_arg.lambda);
          wavefunc[i] = DawsonF(psq)/psq;
        }
        break;
      case WAVEFUNC_TYPE_PAIR1:
        wavefunc[i] = 1.0/(psq+two_body_arg.lambda);
        break;
      case WAVEFUNC_TYPE_PAIR2:
        if (psq < 1e-15) {
          wavefunc[i] = 100.0;
        } else {
          wavefunc[i] = exp(- two_body_arg.lambda*sqrt(psq) )/psq;
        }
        break;
      default:
        ERR.NotImplemented(fname, "WavefuncType not supported.");
    }

      // End mess

  }

}

GeneralizedTwoBody2::~GeneralizedTwoBody2()
{

  const char* fname = "GenerlizedTwoBody2::~GenerlizedTwoBody2()";
  VRB.Func(fname);

  free(wavefunc);
  free(opposite_parity_index);
  free(shift_index);
}




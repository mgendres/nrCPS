#include <iostream>
#include <stdlib.h>
using namespace std;
#include "lattice.h"
#include "arg.h"
#include "random.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

Lattice::Lattice(const LatticeArg &lattice_arg, Random* rng)
{

  rng_p = rng;

  field_type = lattice_arg.field_type;

  field_size = 2*GJP.Vol();
  field = (Float *) malloc(field_size*sizeof(Float)); //Allocate memory for field

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_field, field_size*sizeof(Float)   );
#endif

  Refresh();

}

Lattice::~Lattice()
{

#ifdef USE_GPU
  Cuda::Free( dev_field );
#endif

  free(field);
}


Float* Lattice::Get()
{

#ifdef USE_GPU
  return dev_field;
#else
  return field;
#endif

}

int Lattice::FieldSize()
{
  return field_size;
}

void Lattice::Refresh()
{

  const char* fname = "void Lattice::Refresh()";
  VRB.Func(fname);


  if (field_type == FIELD_TYPE_GAUSS) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = rng_p->Gauss(1.0); // Eventually specify width in Lattice Arg
      field[i+1] = 0.0;
    }
  }

  if (field_type == FIELD_TYPE_COMPLEXGAUSS) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = rng_p->Gauss(1.0); // Eventually specify width in Lattice Arg
      field[i+1] = rng_p->Gauss(1.0); // Eventually specify width in Lattice Arg
    }
  }

  if (field_type == FIELD_TYPE_ZTWO) {
    for(int i=0; i<field_size; i+=2) {
      field[i] = 2.0 * rng_p->Z(2)-1.0; // Eventually specify "p" in Lattice Arg
      field[i+1] = 0.0;
    }
  }

  if (field_type == FIELD_TYPE_ZTHREE) {
    int rand_num;
    for(int i=0; i<field_size; i+=2) {
      rand_num = rng_p->Z(3);
      if (rand_num==0) {
        field[i] = 1.0;
        field[i+1] = 0.0;
      } else if (rand_num==1) {
        field[i] = -0.5;
        field[i+1] = 0.8660254037844388;
      } else {
        field[i] = -0.5;
        field[i+1] = -0.8660254037844388;
      }
    }
  }

#ifdef USE_GPU
  Cuda::MemCopy(dev_field, field, field_size*sizeof(Float), cudaMemcpyHostToDevice);
#endif

}

void Lattice::Refresh(Float shift)
{

  const char* fname = "void Lattice::Refresh(Float)";
  VRB.Func(fname);

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Nope--not at all.");
#endif

  Refresh();
  for ( int i=0; i<GJP.Vol(); ++i)  { field[2*i] += shift; }

}

void Lattice::Refresh(complex<Float> shift)
{

  const char* fname = "void Lattice::Refresh(Float)";
  VRB.Func(fname);

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Nope--not at all.");
#endif

  Float re = real(shift);
  Float im = imag(shift);
  Refresh();
  for ( int i=0; i<GJP.Vol(); ++i)  {
    field[2*i] += re;
    field[2*i+1] += im;
  }

}



void Lattice::Refresh(Float* shift)
{

  const char* fname = "void Lattice::Refresh(Float*)";
  VRB.Func(fname);

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Nope--not at all.");
#endif

  Refresh();
  for ( int i=0; i<2*GJP.Vol(); ++i)  { field[i] += shift[i]; }

}

complex<Float> Lattice::Mean()
{

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Nope--not at all.");
#endif


  Float re = 0.0;
  Float im = 0.0;
  int vol = GJP.Vol();
  for (int i=0; i<vol; ++i) {
    re += field[2*i];
    im += field[2*i+1];
  }
  re /= vol;
  im /= vol;
  return ( complex<Float>(re,im) );
}

Float Lattice::Variance()
{

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Nope--not at all.");
#endif

  complex<Float> mean = Mean();
  Float variance = 0.0;
  int vol = GJP.Vol();
  for (int i=0; i<2*vol; ++i) {
    variance += field[i]*field[i];
  }
  variance /= vol;
  variance -= norm(mean);
  return variance;
}



Float& Lattice::operator[](std::size_t position)
{
  return field[position];
}


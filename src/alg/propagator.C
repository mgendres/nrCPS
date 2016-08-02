#include <iostream>
#include <vector>
#include <math.h>
using namespace std;
#include "dispersion.h"
#include "hamiltonian.h"
#include "propagator.h"
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "utils.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

Propagator::Propagator(const PropagatorArg &propagator_arg)
{

  const char* fname = "Propagator::Propagator(const PropagatorArg &)";

  vol = GJP.Vol();
  b = new Fourier(propagator_arg.n_fermions);

}

Propagator::~Propagator()
{

  delete b;

}

void Propagator::Run(vector<Hamiltonian*> hamiltonians)
{
  const char* fname = "void Propagator::Run(vector<Hamiltonian*>)";
  VRB.Func(fname);

  //---- Apply product ofinteraction matrices, b -> e^{-H}.b
  for (int i=0; i<hamiltonians.size(); ++i) {
    hamiltonians[i]->Evolve(b);
  }

}

void Propagator::Set(Float* psi, int n)
{
  const char* fname = "void Propagator::SetSource(Float*)";
  VRB.Func(fname);

  Float* a = Get(n);

#ifdef USE_GPU
  Cuda::MemCopy(a, psi, 2*vol*sizeof(Float), cudaMemcpyDeviceToDevice);
#else
  for (int i=0; i<2*vol; ++i) {
    a[i] = psi[i];
  }
#endif



}

Float* Propagator::Get(int n)
{

  const char* fname = "Float* Propagator::Get()";
  VRB.Func(fname);
  return (Float*) b->Get(n);

}

Fourier* Propagator::Get() { return b; }


void Propagator::Normalize()
{

  const char* fname = "void Propagator::Normalize()";
  VRB.Func(fname);

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Sorry.");
#endif

//#ifdef USE_GPU
//  Cuda::MemCopy(b, dev_b, 2*vol*sizeof(Float), cudaMemcpyDeviceToHost);
//#endif

  for (int n=0; n<b->BatchSize(); ++n) {

    Float* a = (Float*) b->Get(n);
    Float norm = 0.0;

    for (int i=0; i<2*vol; ++i) { norm += a[i]*a[i]; }
    norm = sqrt(norm);

    for (int i=0; i<2*vol; ++i) { a[i] /= norm; }

  }

//#ifdef USE_GPU
//  Cuda::MemCopy(dev_b, b, 2*vol*sizeof(Float), cudaMemcpyHostToDevice);
//#endif



}

void Propagator::Normalize(Float norm)
{

  const char* fname = "void Propagator::Normalize(Float)";
  VRB.Func(fname);

#ifdef USE_GPU
  ERR.NotImplemented(fname,"Sorry.");
#endif

//#ifdef USE_GPU
//  Cuda::MemCopy(b, dev_b, 2*vol*sizeof(Float), cudaMemcpyDeviceToHost);
//#endif

  for (int n=0; n<b->BatchSize(); ++n) {
    Float* a = (Float*) b->Get(n);
    for (int i=0; i<2*vol; ++i) { a[i] /= norm; }
  }

//#ifdef USE_GPU
//  Cuda::MemCopy(dev_b, b, 2*vol*sizeof(Float), cudaMemcpyHostToDevice);
//#endif

}


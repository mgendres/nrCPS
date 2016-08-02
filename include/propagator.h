#include "config.h"
#include <vector>
#include "hamiltonian.h"
#include "arg.h"
#include "fourier.h"

#ifndef INCLUDED_PROP
#define INCLUDED_PROP

class Propagator
{
  private:
    Fourier* b;
    int vol;

    Propagator& operator=(const Propagator&);
    Propagator(const Propagator&);

  public:
    explicit Propagator(const PropagatorArg &);
    ~Propagator();
    void Run( vector<Hamiltonian*> );
    void Set(Float*, int);
    Float* Get(int);
    Fourier* Get();
    void Normalize();
    void Normalize(Float);
};

#endif

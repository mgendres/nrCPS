#include <complex>
#include "config.h"
#include "arg.h"
#include "random.h"

#ifndef INCLUDED_LATTICE
#define INCLUDED_LATTICE
class Lattice
{
  private:
    Random* rng_p;
    enum FieldType field_type;
    int field_size;
    Float* field;

#ifdef USE_GPU
    Float* dev_field;
#endif

    Lattice& operator=(const Lattice&);
    Lattice(const Lattice&);

  public:
    explicit Lattice(const LatticeArg &, Random*);
    ~Lattice();
    void Refresh();
    void Refresh(Float);
    void Refresh(complex<Float>);
    void Refresh(Float*);
    Float* Get();
    int FieldSize();
    complex<Float> Mean();
    Float Variance();
    Float& operator[](std::size_t);
};
#endif

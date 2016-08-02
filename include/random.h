#include "arg.h"
#include "uniform_deviates.h"

#ifndef INCLUDED_RANDOM
#define INCLUDED_RANDOM
class Random
// Random number generator base class
{
  private:
    UniformDeviate* r;
    RandomType random_type; // This can be made obsolete...

    Random& operator=(const Random&);
    Random(const Random&);

  public:
    explicit Random(const RandomArg &);
    ~Random();
    Float Uniform();
    Float Gauss(Float);
    Float Normal(Float, Float);
    Float LogNormal(Float, Float);
    int Z(int);
    void ReadState(const char*);
    void WriteState(const char*);
};

#endif

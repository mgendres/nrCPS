#include "enum.h"

#ifndef INCLUDED_DISPERSION
#define INCLUDED_DISPERSION
class Dispersion
{
  private:
    Float* dispersion;

    Dispersion& operator=(const Dispersion&);
    Dispersion(const Dispersion&);

  public:
    explicit Dispersion(DispersionType, Float, CutoffType);
    ~Dispersion();
    Float Get(int, int, int); // Returns energy associated with momenta specifies by three integers
    Float Get(int); // Returns energy associated with momenta specifies by collective index
};

#endif

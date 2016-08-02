#include <complex>
#include "one_body.h"
#include "two_body.h"
#include "enum.h"

#ifndef INCLUDED_SLATER_DET2
#define INCLUDED_SLATER_DET2

class SlaterDet2
{
  private:
    complex<Float>* phi;  // The NxN matrix:  <phi|i_1, i_2>
    complex<Float>* psi1; // The NxN matrix:  <f_1| i_1>
    complex<Float>* psi2; // The NxN matrix:  <f_2| i_2>
    int N;                 // Total number of fermions
    DetType det_type;
    void Initialize(int);

    SlaterDet2& operator=(const SlaterDet2&);
    SlaterDet2(const SlaterDet2&);

  public:
    explicit SlaterDet2(Float**, Float**, OneBody**, OneBody**, TwoBody*, int);
    ~SlaterDet2();
    complex<Float> Run(int, int, int);
    complex<Float> SlaterMat(int, int);
    void SetDetType(DetType);
};

#endif

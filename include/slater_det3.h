#include <complex>
#include "one_body.h"
#include "two_body.h"
#include "enum.h"

#ifdef USE_GPU
#endif

#ifndef INCLUDED_SLATER_DET3
#define INCLUDED_SLATER_DET3

class SlaterDet3
{
  private:
    complex<Float>* phi;  // The NxN matrix:  <phi|i_1, i_2>
    complex<Float>* psi1; // The NxN matrix:  <f_1| i_1>
    complex<Float>* psi2; // The NxN matrix:  <f_2| i_2>
    int N;                 // Total number of fermions
    DetType det_type;

#ifdef USE_GPU
    Float** props1;
    Float** props2;
    Float** one_body1;
    Float** one_body2;
#else
    Float** props1;
    Float** props2;
    OneBody** one_body1;
    OneBody** one_body2;
#endif

    SlaterDet3& operator=(const SlaterDet3&);
    SlaterDet3(const SlaterDet3&);

  public:
    explicit SlaterDet3(Float**, Float**, OneBody**, OneBody**, int);
    ~SlaterDet3();
    complex<Float> Run(int, int, int);
    complex<Float> SlaterMat(int, int);
    void Initialize(TwoBody*);
    void SetDetType(DetType);
};

#endif

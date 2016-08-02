#include <complex>
#include "one_body.h"
#include "two_body.h"
#include "enum.h"

#ifndef INCLUDED_SLATER_DET
#define INCLUDED_SLATER_DET

class SlaterDet
{
  private:
    complex<Float>* slater_mat;
    int mat_size;
    DetType det_type;
    void Initialize(int);

    SlaterDet& operator=(const SlaterDet&);
    SlaterDet(const SlaterDet&);

  public:
    explicit SlaterDet(Float**, OneBody**, int);
    explicit SlaterDet(Float**, Float**, TwoBody*, int);
    ~SlaterDet();
    complex<Float> Run(int);
    complex<Float> SlaterMat(int, int);
    void SetDetType(DetType);
};

#endif

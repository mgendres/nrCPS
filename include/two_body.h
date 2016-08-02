#include <complex>
#include "config.h"
#include "arg.h"
#include "vector3.h"
#include "oh_group.h"

#ifndef INCLUDED_TWO_BODY
#define INCLUDED_TWO_BODY

class TwoBody
{
  private:
    Float* wavefunc;
    int* opposite_parity_index;
    int vol;
    Vector3 dims;
    Vector3 bcs;

#ifdef USE_GPU
    Float* dev_wavefunc;
    int* dev_opposite_parity_index;
    Float* psum;
    Float* dev_psum;
    int dim;
#endif

    TwoBody& operator=(const TwoBody&);
    TwoBody(const TwoBody&);

  public:
    explicit TwoBody(const TwoBodyArg &);
    ~TwoBody();
    void Deform(string); // Set source based on parameters specified in a file
    complex<Float> Run(Float*, Float*);
    void Run(Float**, Float**, complex<Float>*, int);
};



class GeneralizedTwoBody
{
  private:
    Float* half_site_wavefunc; // 2^d times the size of wave func, d = space-dim
    int* opposite_parity_index;
    int* shift_indices;
    int vol;
    Vector3 dims;
    Vector3 bcs;

    OhGroup oh_group;

    GeneralizedTwoBody& operator=(const GeneralizedTwoBody&);
    GeneralizedTwoBody(const GeneralizedTwoBody&);

  public:
    explicit GeneralizedTwoBody(const GeneralizedTwoBodyArg &);
    ~GeneralizedTwoBody();
    complex<Float> Run(Float*, Float*, OctahedralRep);
    complex<Float> Run2(Float*, Float*, OctahedralRep);

};



class GeneralizedTwoBody2
{
  private:
    Float* wavefunc;
    int* opposite_parity_index;
    int* shift_index;
    int vol;
    Vector3 dims;
    Vector3 bcs;
    Vector3 shift;

    GeneralizedTwoBody2& operator=(const GeneralizedTwoBody2&);
    GeneralizedTwoBody2(const GeneralizedTwoBody2&);

  public:
    explicit GeneralizedTwoBody2(const GeneralizedTwoBodyArg &);
    ~GeneralizedTwoBody2();
    complex<Float> Run(Float*, Float*);

};




#endif

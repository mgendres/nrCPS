#include <complex>
#include "config.h"
#include "arg.h"
#include "vector3.h"

#ifndef INCLUDED_ONE_BODY
#define INCLUDED_ONE_BODY

class OneBody
{
  private:
    SourceType source_type;
    Float lambda1;
    Float lambda2;
    Float lambda3;
    Float* psi; // Source vector
    int vol;

#ifdef USE_GPU
    Float* dev_psi;
    Float* psum;
    Float* dev_psum;
    int dim;
#endif

    OneBody& operator=(const OneBody&);
    OneBody(const OneBody&);

  public:
    explicit OneBody(const OneBodyArg &);
    ~OneBody();
    void Set(string); // Set source based on parameters specified in a file
    void Set(); // Set source based on what is stored in psi and specified source_type.
    void Set(int, int, int); // Set source to basis vector labeled by ( coord1, coord2, coord3)
    void Set(int, int, int, int); // Set source to basis vector labeled by ( coord1, coord2, coord3)
    void Set(Float*); // Set source to external vector
    void Set(const Vector3&);
    void Set(const Vector3&, OctahedralRep);
    complex<Float> Project(Float*); // perform projection of argument vertor onto source vector
    Float* Get(); // Get the address to the source vector
};

#endif

#include "config.h"
#include "enum.h"
#include <vector3.h>


#ifndef INCLUDED_OH_GROUP
#define INCLUDED_OH_GROUP

class OhGroup
{
  private:
    // These are just useful tables
    static const int chi[48][10];
    static const int dim[10];
    static const int subshell_dimensions[8];

    // These are stuff computed in the constructor
    Vector3 vectors[48];
    Float weights[48][10];
    int reduced;
    int subshell_texture;
    int subshell_dimension;
    int* opposite_parity_index;

    // Some private class functions
    Vector3 GetDistinctVector(const Vector3&);
    SubshellTexture GetSubshellTexture(const Vector3&);

    OhGroup& operator=(const OhGroup&);
    OhGroup(const OhGroup&);

  public:
    explicit OhGroup(const Vector3&, int);
    ~OhGroup();

    int GetSubshellDimension();
    int GetSubshellTexture();

    Vector3 GetVector(int);
    Vector3 GetVector(int, Parity);
    int GetOppositeParityIndex(int);
    Float GetWeight(OctahedralRep, int);

};



class OhProjector
{

  private:
    int subspace_dimension;
    int u_dim;
    Float* u;

  public:
    explicit OhProjector(const Vector3&, OctahedralRep);
    explicit OhProjector(const Vector3&, const Vector3&, OctahedralRep);
    explicit OhProjector(const Vector3&, const Vector3&, const Vector3&, OctahedralRep);
    ~OhProjector();
  
    int GetProjectorDimension();
    int GetSubspaceDimension();
    Float* GetProjectorEigenvector(int);

};


#endif

#include <string>
#include "arg.h"
#include "vector3.h"

// Header file for the verbose class. An object of this class
// called VRB should be created at the highest scope (outside 
// main). The header file declares VRB as external.

#ifndef INCLUDED_GLOBAL_JOB_PARAMETER
#define INCLUDED_GLOBAL_JOB_PARAMETER
class GlobalJobParameter
{
  private:
    int t_sites;
    int x_sites;
    int y_sites;
    int z_sites;
    BoundaryType x_boundary_type;
    BoundaryType y_boundary_type;
    BoundaryType z_boundary_type;
    Float cutoff; 

  public:
    explicit GlobalJobParameter();
    ~GlobalJobParameter();
    void Initialize(const DoArg &);
    int Tsites();
    int Xsites();
    int Ysites();
    int Zsites();
    int Vol();
    Vector3 GetDimensions();
    BoundaryType Xbc();
    BoundaryType Ybc();
    BoundaryType Zbc();
    bool APBCQuery();
    Vector3 GetBoundaryConditions();
    void GetBcShift(Float*); // used to find shift in momenta due to bc p_n = 2 Pi (n + shift)/ L
    Float Cutoff();
    Float CutoffSq();
};

extern GlobalJobParameter GJP;
#endif

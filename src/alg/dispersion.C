#include <iostream>
#include <limits>
#include <math.h>
#include <stdlib.h>
using namespace std;
#include "dispersion.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"

Dispersion::Dispersion(DispersionType dispersion_type, Float mass, CutoffType cutoff_type)
{

  const char* fname = "Dispersion::Dispersion(DispersionType)";

  dispersion = (Float *) malloc(GJP.Vol()*sizeof(Float));  

  Float bc[3];
  GJP.GetBcShift(bc);

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  Float sx;
  Float sy;
  Float sz;

  Float psq;

  int z;
  int y;
  int x;
  Vector3 vec;

  for(int i=0; i<GJP.Vol(); ++i) {

    vec = IndexToVector3(i);
    x = vec.GetX();
    y = vec.GetY();
    z = vec.GetZ();

    if (dispersion_type==DISPERSION_TYPE_STANDARD) {
      // \Delta(p) = 2 \sum_j \sin^2(p_j/2)
      sx = sin( (x+bc[0]) * PI / x_sites );
      sy = sin( (y+bc[1]) * PI / y_sites );
      sz = sin( (z+bc[2]) * PI / z_sites );
      dispersion[i] = 2.0 * ( sx*sx + sy*sy + sz*sz ); 
      dispersion[i] /= mass; 
    }
  
    if ( (dispersion_type==DISPERSION_TYPE_QUADRATIC)||(dispersion_type==DISPERSION_TYPE_PERFECT) ) {

      if (x >= x_sites/2) { x = x_sites - x - GJP.Xbc(); }
      if (y >= y_sites/2) { y = y_sites - y - GJP.Ybc(); }
      if (z >= z_sites/2) { z = z_sites - z - GJP.Zbc(); }

      sx = (x + bc[0]) / (Float)x_sites;
      sy = (y + bc[1]) / (Float)y_sites;
      sz = (z + bc[2]) / (Float)z_sites;
      psq = (sx*sx + sy*sy + sz*sz)*4.0*PISQ;

      //---- Cutoff, beyond which we take psq = infinity
      if ( (cutoff_type==CUTOFF_TYPE_HARD)&&(psq > GJP.CutoffSq() ) ) {
        dispersion[i] = numeric_limits<Float>::infinity();
      } else {
        if (dispersion_type==DISPERSION_TYPE_QUADRATIC) { dispersion[i] = psq/(2.0*mass); }
        if (dispersion_type==DISPERSION_TYPE_PERFECT) { dispersion[i] = ( exp( psq/(2.0*mass) ) -1.0 ); }
      }

    }

  }

}

Dispersion::~Dispersion()
{
  free(dispersion);
}


Float Dispersion::Get(int n1, int n2, int n3)
{
  const char* fname = "Float Dispersion::Get(int, int, int)";
  return dispersion[n3+n2*GJP.Zsites()+n1*GJP.Ysites()*GJP.Zsites()];
}

Float Dispersion::Get(int i)
{
  const char* fname = "Float Dispersion::Get(int)";
  return dispersion[i];
}

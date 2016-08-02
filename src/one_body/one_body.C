#include <iostream>
#include <fstream>
#include <string>
#include <math.h>
#include <stdlib.h>

using namespace std;
#include "arg.h"
#include "one_body.h"
#include "verbose.h"
#include "special_functions.h"
#include "error.h"
#include "global_job_parameter.h"
#include "constants.h"
#include "oh_group.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

OneBody::OneBody(const OneBodyArg &one_body_arg)
{

  const char* fname = "void OneBody::OneBody(const OneBodyArg &)";
  VRB.Func(fname);

  source_type = one_body_arg.source_type;
  lambda1 = one_body_arg.lambda1; 
  lambda2 = one_body_arg.lambda2; 
  lambda3 = one_body_arg.lambda3; 
  vol = GJP.Vol();
  psi = (Float *) malloc(2*vol*sizeof(Float));  

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_psi, 2*vol*sizeof(Float) );

  psum = (Float *) malloc(2*BLOCKS*sizeof(Float));
  Cuda::Malloc( (void**)&dev_psum, 2*BLOCKS*sizeof(Float) );
#endif

}

OneBody::~OneBody()
{
#ifdef USE_GPU
  free(psum);
  Cuda::Free( dev_psum );
  Cuda::Free( dev_psi );
#endif
  free(psi);
}

void OneBody::Set()
{

  const char* fname = "void OneBody::Set()";
  VRB.Func(fname);

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();
  int vol = GJP.Vol();

  if (source_type==SOURCE_TYPE_MOM) {

    //---- No need to do anything!

  }

  if (source_type==SOURCE_TYPE_SHO) {

    if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"Must use PBCs when a potential is present."); }

    //---- This is the SHO state in momentum space

    int j;
    int n1_max = 10;
    int n2_max = 10;
    int n3_max = 10;
  
    if ( n1_max > x_sites ) {
      n1_max = x_sites;
      VRB.Warn(fname, "Reset maximum SHO level nx to %d", n1_max);
    }
    if ( n2_max > y_sites ) {
      n2_max = y_sites;
      VRB.Warn(fname, "Reset maximum SHO level ny to %d", n2_max);
    }
    if ( n3_max > z_sites ) {
      n3_max = z_sites;
      VRB.Warn(fname, "Reset maximum SHO level nz to %d", n2_max);
    }

    VRB.Warn(fname, "Maximum SHO levels (nx,ny,nz) are (%d,%d,%d)",n1_max,n2_max,n3_max);

    //---- Copy coeffs from psi to coeffs
    Float coeffs[n1_max][n2_max][n3_max][2];
    for (int n1=0; n1<n1_max; ++n1) {
      for (int n2=0; n2<n2_max; ++n2) {
        for (int n3=0; n3<n3_max; ++n3) {
          j = 2*Vector3ToIndex( Vector3(n1,n2,n3) );
          coeffs[n1][n2][n3][0] = psi[j];
          coeffs[n1][n2][n3][1] = psi[j+1];
        }
      }
    }
 
    //---- Now lets set the source
    Float dx, dy, dz;
    int x, y, z;
    Float h1[n1_max];
    Float h2[n2_max];
    Float h3[n3_max];
    Vector3 vec;

    for(int i=0; i<2*vol;i+=2) {

      //---- Determine coordinates
      vec = IndexToVector3( i/2 );
      x = vec.GetX();
      y = vec.GetY();
      z = vec.GetZ();

      //---- Momentum space wve function is centered around p=0 (PBCS)
      if (x>=x_sites/2) { x -= x_sites; }
      if (y>=y_sites/2) { y -= y_sites; }
      if (z>=z_sites/2) { z -= z_sites; }

      //---- For each set of integers, compute p_j = 2 pi k_j / L_j and multiply by (L_0)_j
      dx = TWOPI * x * lambda1 / (Float) x_sites;
      dy = TWOPI * y * lambda2 / (Float) y_sites;
      dz = TWOPI * z * lambda2 / (Float) z_sites;

      //---- Compute Hermite polynomials
      HermiteH(h1, n1_max, dx);
      HermiteH(h2, n2_max, dy);
      HermiteH(h3, n3_max, dz);

      //---- Compute wave function as a linear combination of SHO basis wave functions
      Float tmp;
      psi[i] = 0.0;
      psi[i+1] = 0.0; 
      for (int n1=0; n1<n1_max; ++n1) {
        for (int n2=0; n2<n2_max; ++n2) {
          for (int n3=0; n3<n3_max; ++n3) {
            tmp = h1[n1]*h2[n2]*h3[n3];
            psi[i] += coeffs[n1][n2][n3][0] * tmp;
            psi[i+1] += coeffs[n1][n2][n3][1] * tmp; 
          }
        }
      }
      
      tmp = exp( -(dx*dx+dy*dy+dz*dz)/2.0 );
      tmp *= sqrt(lambda1*lambda2*lambda3);

      //---- Finally, scale result by (2pi/L)^3
      tmp *= sqrt( (TWOPI*TWOPI*TWOPI)/vol );

      psi[i] *= tmp;
      psi[i+1] *= tmp; 

    }


  }

#ifdef USE_GPU
  Cuda::MemCopy(dev_psi, psi, 2*vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif

}

void OneBody::Set(string vml_filename)
{

  const char* fname = "void OneBody::Set(string)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; ++i) { psi[i]=0.0; }

  ifstream file;
  file.open(vml_filename.c_str());
  if (!file) { ERR.FileR(fname,vml_filename.c_str()); }

  int coord1;
  int coord2;
  int coord3;

  int j;
  //---- Need to add error checking, etc....
  while (true) {

    file >> coord1;
    file >> coord2;
    file >> coord3;

    j = 2*Vector3ToIndex( Vector3(coord1, coord2, coord3) );

    if ( j < 2*GJP.Vol() ) {
      file >> psi[j]; 
      file >> psi[j+1];
    } else {
      ERR.General(fname, "Source quantum numbers excede lattice vector size.");
    }

    if ( file.eof() ) break; // Should be moved two lines up?--if EOF, then accessing random memory.
  }

  file.close();

  Set(); 

}

void OneBody::Set(int coord1, int coord2, int coord3)
{

  const char* fname = "void OneBody::Set(int, int, int)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; ++i) { psi[i] = 0.0; }

  int j = Vector3ToIndex( Vector3(coord1, coord2, coord3) );
  psi[2*j] = 1.0; 

  Set(); 

}

void OneBody::Set(int coord1, int coord2, int coord3, int parity)
{

  const char* fname = "void OneBody::Set(int, int, int)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; ++i) { psi[i] = 0.0; }

  int j = Vector3ToIndex( Vector3(coord1, coord2, coord3) );
  psi[2*j] = 1.0; 

  j = Vector3ToIndex( -1 * Vector3(coord1, coord2, coord3) );

  if (parity) {
     psi[2*j] = 1.0;
  } else {
     psi[2*j] = -1.0;
  } 

  Set(); 

}



void OneBody::Set(Float* source) {

  const char* fname = "void OneBody::Set(Float*)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; ++i) { psi[i] = source[i]; }

#ifdef USE_GPU
  Cuda::MemCopy(dev_psi, psi, 2*vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif

}

void OneBody::Set(const Vector3& vec)
{

  const char* fname = "void OneBody::Set(const Vector3&)";
  VRB.Func(fname);

  for (int i=0; i<2*vol; ++i) { psi[i] = 0.0; }

  int j = Vector3ToIndex(  vec );
  psi[2*j] = 1.0;

  Set(); 

}



void OneBody::Set(const Vector3& vec, OctahedralRep rep)
{

  const char* fname = "void OneBody::Set(const Vector3&, OctahedralRep)";
  VRB.Func(fname);

  int j;
  OhGroup oh_group(vec, 1);

  for (int i=0; i<2*vol; ++i) { psi[i] = 0.0; }

  // Assumes vec is a momentum vector!!!
  for (int l=0; l<oh_group.GetSubshellDimension(); ++l) {
    j = Vector3ToIndex(  oh_group.GetVector(l) );
    psi[2*j] += oh_group.GetWeight(rep, l) ;
  }

  Set(); 

}


Float* OneBody::Get()
{

#ifdef USE_GPU
  return dev_psi;
#else
  return psi;
#endif

}

#include <iostream>
#include <math.h>
#include <stdlib.h>

using namespace std;
#include "lattice.h"
#include "arg.h"
#include "verbose.h"
#include "error.h"
#include "global_job_parameter.h"
#include "hamiltonian.h"
#include "dispersion.h"
#include "constants.h"
#include "fourier.h"

#ifdef USE_GPU
#include "cuda_utils.h"
#endif

Hamiltonian::Hamiltonian()
{

  const char* fname = "void Hamiltonian::Hamiltonian()";
  VRB.Func(fname);

}

Hamiltonian::~Hamiltonian()
{

  const char* fname = "void Hamiltonian::~Hamiltonian()";
  VRB.Func(fname);

}

Interactions::Interactions(  vector<Hamiltonian*>  interaction_list)
{

  const char* fname = "Interactions::Interactions(  vector<Hamiltonian*> )";
  VRB.Func(fname);

  vol = GJP.Vol();

#ifdef USE_GPU
  Cuda::Malloc( (void**)&a, 2*vol*sizeof(Float) );
  Cuda::Malloc( (void**)&c, 2*vol*sizeof(Float) );
#else
  a = (Float *) malloc(2*vol*sizeof(Float));  
  c = (Float *) malloc(2*vol*sizeof(Float));  
#endif

  interactions = interaction_list;

}

void Interactions::Initialize() {

  for (int n=0; n<interactions.size(); ++n) { interactions[n]->Initialize(); }

}



Interactions::~Interactions()
{

#ifdef USE_GPU
  Cuda::Free(a);
  Cuda::Free(c);
#else
  free(a);
  free(c);
#endif

}

Interaction::Interaction(Lattice* lat, const InteractionArg &interaction_arg, Float dt)
: lattice(lat)
{

  const char* fname = "Interaction::Interaction(Lattice*, const InteractionArg &, Float)";
  VRB.Func(fname);

  //---- Compute the interactions in momentum space
  int vol = GJP.Vol();
  interaction = (Float *) malloc(vol*sizeof(Float));  


  Float mass = interaction_arg.mass;

  Dispersion dispersion(DISPERSION_TYPE_PERFECT, mass, CUTOFF_TYPE_NONE);
  Float xi;

  Float psq;

  Float PSQ = GJP.CutoffSq();
  Float XI = exp( PSQ/(2.0*mass));

  for(int i=0; i<vol; ++i) {

    xi = 1.0 + dispersion.Get(i);
    if (xi > XI) { xi = XI; }

    switch (interaction_arg.interaction_type) {
      case INTERACTION_TYPE_NONE:
        ERR.NotImplemented(fname,"Interaction type INTERACTION_TYPE_NONE.");
        break;
      case INTERACTION_TYPE_ONEMINUSXIINVSQ:
        psq = 1.0 - 1.0 /(xi*xi); 
        break;
      case INTERACTION_TYPE_XISQMINUSONE:
        psq = xi*xi - 1.0; 
        break;
      default:
        ERR.NotImplemented(fname,"Unrecognized interaction_type.");
    }

    psq *= mass;

    //---- Evaluate the interaction as a Taylor series in p^2: O(p) = sum_n C_n p^(2n)
    interaction[i] = 0.0;
    for (int j=0; j<interaction_arg.num_couplings; ++j) {
      interaction[i] += interaction_arg.couplings[j]*pow(psq,j);
    }

    //---- Make sure O(p) is non-negative before taking the square root!
    if (interaction[i]<0.0) {
      ERR.General(fname, "Interaction is less the zero; cannot take square root.");
    }

    //---- Normalize operator, etc..
    interaction[i] *= 2.0 * TWOPI / mass;  // Just a convention
    interaction[i] = sqrt(interaction[i]);
    interaction[i] /= vol; // Division by vol because FFT is not normalized
//    cout << i << " " << interaction[i] << "\n";
    interaction[i] *= dt; // Controls overall sign and step size of the interaction

  }

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_interaction, vol*sizeof(Float) );
  Cuda::MemCopy(dev_interaction, interaction, vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif

  field = new Fourier(1);  

}

Interaction::~Interaction()
{

  const char* fname = "void Interaction::~Interaction()";
  VRB.Func(fname);

  delete field;

#ifdef USE_GPU
  Cuda::Free( dev_interaction );
#endif

  free(interaction);

}

Potential::Potential(const PotentialArg &potential_arg, Float dt)
{

  const char* fname = "Potential::Potential(const PotentialArg &, Float)";
  VRB.Func(fname);

  if ( GJP.APBCQuery() ) { ERR.NotImplemented(fname,"Must use PBCs when a potential is present."); }

  int vol = GJP.Vol();
  potential = (Float *) malloc(vol*sizeof(Float));  

  //---- Compute the potential
  int x;
  int y;
  int z;
  int Y; // Y = y+x*Y

  int x_sites = GJP.Xsites();
  int y_sites = GJP.Ysites();
  int z_sites = GJP.Zsites();

  Float v; 

  for(int i=0; i<vol; ++i) {

    z = i%z_sites;
    Y = i/z_sites;
    y = Y%y_sites;
    x = Y/y_sites;

    if (x>=x_sites/2) { x -= x_sites; }
    if (y>=y_sites/2) { y -= y_sites; }
    if (z>=z_sites/2) { z -= z_sites; }

    //----- Potential form
    switch (potential_arg.potential_form) {
      case POTENTIAL_FORM_NONE:
        v = 0.0;
        break;
      case POTENTIAL_FORM_HARMONIC:
        v =  potential_arg.spring_constant1*x*x;
        v += potential_arg.spring_constant2*y*y;
        v += potential_arg.spring_constant3*z*z;
        v /= 2.0;
        break;
      case POTENTIAL_FORM_COULOMB:
        ERR.NotImplemented(fname, "Potential form POTENTIAL_FORM_COULOMB not supported.");
        v =  x*x/potential_arg.spring_constant1;
        v += y*y/potential_arg.spring_constant2;
        v += z*z/potential_arg.spring_constant3;
        v = 1.0/sqrt(v);
        break;
      default:
        ERR.NotImplemented(fname, "Potential form not supported.");
    }

    //---- Potential discretization
    switch (potential_arg.potential_type) {
      case POTENTIAL_TYPE_LIN:
        potential[i] = 1.0 - dt*v;
        potential[i] /= vol;
        break;
      case POTENTIAL_TYPE_EXP:
        potential[i] = exp(-dt*v);
        potential[i] /= vol;
        break;
      default:
        ERR.NotImplemented(fname, "Potential type not supported.");
    }

  }

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_potential, vol*sizeof(Float) );
  Cuda::MemCopy(dev_potential, potential, vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif



}

Potential::~Potential()
{

  const char* fname = "void SplitPotential::~Potential()";
  VRB.Func(fname);

#ifdef USE_GPU
  Cuda::Free( dev_potential );
#endif

  free(potential);

}

Kinetic::Kinetic(const KineticArg &kinetic_arg, Float dt)
{

  const char* fname = "Kinetic::Kinetic(const KineticArg &, Float)";
  VRB.Func(fname);

  int vol = GJP.Vol();
  xi = (Float *) malloc(vol*sizeof(Float));

  Dispersion dispersion(kinetic_arg.dispersion_type, kinetic_arg.mass/dt, CUTOFF_TYPE_HARD);
  for (int i=0; i<vol; ++i) {
    xi[i] = 1.0 + dispersion.Get(i);
  }

#ifdef USE_GPU
  Cuda::Malloc( (void**)&dev_xi, vol*sizeof(Float) );
  Cuda::MemCopy(dev_xi, xi, vol*sizeof(Float), cudaMemcpyHostToDevice);
#endif



}

Kinetic::~Kinetic()
{

  const char* fname = "void Kinetic::~Kinetic()";
  VRB.Func(fname);

#ifdef USE_GPU
  Cuda::Free( dev_xi );
#endif

  free(xi);

}


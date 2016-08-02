#include "arg.h"
#include "lattice.h"
#include "fourier.h"

#ifndef INCLUDED_HAMILTONIAN
#define INCLUDED_HAMILTONIAN

// Abstract class from which HamiltonianYukawa and Potential are derived
class Hamiltonian
{
  private:

    Hamiltonian& operator=(const Hamiltonian&);
    Hamiltonian(const Hamiltonian&);

  public:
    explicit Hamiltonian();
    virtual ~Hamiltonian();
    virtual void Evolve(Fourier*) = 0;
    virtual void Initialize() = 0;
};

// Classes derived from Hamiltonian 

class Interactions : public Hamiltonian
{
  private:
    int vol;
    Float *a;
    Float *c;
    vector<Hamiltonian*> interactions;

    Interactions& operator=(const Interactions&);
    Interactions(const Interactions&);

  public:
    explicit Interactions( vector<Hamiltonian*> );
    ~Interactions();
    void Evolve(Fourier*);
    void Initialize();
};

class Interaction : public Hamiltonian
{
  private:
    Lattice* lattice;
    Fourier* field;
    Float* interaction;
    Float* dev_interaction;

    Interaction& operator=(const Interaction&);
    Interaction(const Interaction&);

  public:
    explicit Interaction(Lattice*, const InteractionArg &, Float);
    ~Interaction();
    void Evolve(Fourier*);
    void Initialize();
};

class Potential : public Hamiltonian
{
  private:
    Float* potential;
    Float* dev_potential;
    PotentialForm potential_form;
    PotentialType potential_type;

    Potential& operator=(const Potential&);
    Potential(const Potential&);

  public:
    explicit Potential(const PotentialArg &, Float);
    ~Potential();
    void Evolve(Fourier*);
    void Initialize();
};

class Kinetic : public Hamiltonian
{
  private:
    Float* xi;
    Float* dev_xi;

    Kinetic& operator=(const Kinetic&);
    Kinetic(const Kinetic&);

  public:
    explicit Kinetic(const KineticArg &, Float);
    ~Kinetic();
    void Evolve(Fourier*);
    void Initialize();
};



#endif

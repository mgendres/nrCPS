#include <string>
#include <vector>
#include <fstream>
#include "enum.h"
#include "vector3.h"

#ifndef INCLUDED_ARG
#define INCLUDED_ARG

// Abstract class from which all Arg classes are derived
class Arg
{
  private:
    std::ifstream file;
    std::string input;
  public:
    explicit Arg();
    virtual ~Arg();
    virtual void Decode(std::string)=0; // Function that imports parameters from a file
    void OpenFile(std::string);
    void CloseFile();

    template <typename T>
    int ReadLine(T &);
    int ReadLine(bool &);
    int ReadLine(RandomType &); // Can't use a template for these...
    int ReadLine(SeedType &);
    int ReadLine(FieldType &);
    int ReadLine(SourceType &);
    int ReadLine(WavefuncType &);
    int ReadLine(BoundaryType &);
    int ReadLine(DispersionType &);
    int ReadLine(CutoffType &);
    int ReadLine(PotentialType &);
    int ReadLine(PotentialForm &);
    int ReadLine(DetType &);
    int ReadLine(InteractionType &);
    int ReadLine(IoType &);
    int ReadLine(Vector3 &);
    int ReadLine(OctahedralRep &);

};

class DoArg : public Arg
{
  public:
    explicit DoArg();
    explicit DoArg(std::string);
    ~DoArg();
    void Decode(std::string);
    int t_sites; // number of sites in the t direction
    int x_sites; // number of sites in the x direction
    int y_sites; // number of sites in the y direction
    int z_sites; // number of sites in the z direction
    BoundaryType x_boundary_type; // boundary conditions in the x direction for the fermion
    BoundaryType y_boundary_type; // boundary conditions in the y direction for the fermion
    BoundaryType z_boundary_type; // boundary conditions in the z direction for the fermion
    Float cutoff; // Hard momentum cutoff in units of pi
};

class LatticeArg : public Arg
{
  public:
    explicit LatticeArg();
    explicit LatticeArg(std::string);
    ~LatticeArg();
    void Decode(std::string);
    enum FieldType field_type; // type of two-body auxiliary field
};

class RandomArg : public Arg
{
  public:
    explicit RandomArg();
    explicit RandomArg(std::string);
    ~RandomArg();
    void Decode(std::string);
    enum RandomType random_type; // type of random number generator
    enum SeedType seed_type; // type of seed for the random number generator
    long seed; // input value for the seed
    std::string file_stem; // file stem for random state 
};

class PropagatorArg : public Arg
{
  public:
    explicit PropagatorArg();
    explicit PropagatorArg(std::string);
    ~PropagatorArg();
    void Decode(std::string);
    int n_fermions; // number of fermions
    std::string file_stem; // file name stem
};

class OneBodyArg : public Arg
{
  public:
    explicit OneBodyArg();
    explicit OneBodyArg(std::string);
    ~OneBodyArg();
    void Decode(std::string);
    SourceType source_type; // one-particle basis state type for the source
    Float lambda1; // wavefunction parameter
    Float lambda2; // wavefunction parameter
    Float lambda3; // wavefunction parameter
};

class TwoBodyArg : public Arg
{
  public:
    explicit TwoBodyArg();
    explicit TwoBodyArg(std::string);
    ~TwoBodyArg();
    void Decode(std::string);
    WavefuncType wavefunc_type; // two-body wave function type
    DispersionType dispersion_type1;
    DispersionType dispersion_type2;
    Float mass1; // wave function mass parameter for first particle
    Float mass2; // wave function mass parameter for second particle
    Float lambda; // additional wave function parameter 
    std::string file_stem; // file name stem
};

class GeneralizedTwoBodyArg : public Arg
{
  public:
    explicit GeneralizedTwoBodyArg();
    explicit GeneralizedTwoBodyArg(std::string);
    ~GeneralizedTwoBodyArg();
    void Decode(std::string);
    WavefuncType wavefunc_type; // two-body wave function type
//    OctahedralRep rep;
    Vector3 distinct_vec; // distinct vector
    Float lambda; // wave function parameter 
    std::string file_stem; // file name stem
};


class VerboseArg : public Arg
{
  public:
    explicit VerboseArg();
    explicit VerboseArg(std::string);
    ~VerboseArg();
    void Decode(std::string);
    bool func_level;   // function level switch
    bool warn_level;   // warning level switch
    bool result_level; // result level switch
    bool flow_level;   // flow level switch
    bool debug_level;  // debug level switch
};

class EvoArg : public Arg
{
  public:
    explicit EvoArg();
    explicit EvoArg(std::string);
    ~EvoArg();
    void Decode(std::string);
    int start;          // configuration start number
    int unload_period;  // configuration unload period
    int configurations; // total number of configurations to be generated
};

class MomentaArg : public Arg
{
  public:
    explicit MomentaArg();
    explicit MomentaArg(std::string);
    ~MomentaArg();
    void Decode(std::string);
    DispersionType dispersion_type; // dispersion relation used to compute fermi-energy
    Float fermi_energy; // fermi energy
    Float mass; // fermion mass used to compute momenta below fermi energy
};

class InteractionArg : public Arg
{
  public:
    explicit InteractionArg();
    explicit InteractionArg(std::string);
    ~InteractionArg();
    void Decode(std::string);
    Float mass;
    InteractionType interaction_type;
    std::vector<Float> couplings;
    int num_couplings;
};

class PotentialArg : public Arg
{
  public:
    explicit PotentialArg();
    explicit PotentialArg(std::string);
    ~PotentialArg();
    void Decode(std::string);
    PotentialForm potential_form; // Potential functional form (none, harmonic, coulomb)
    PotentialType potential_type; // potential discretization type
    Float spring_constant1;       // harmonic potential spring constant in x-direction
    Float spring_constant2;       // harmonic potential spring constant in y-direction
    Float spring_constant3;       // harmonic potential spring constant in z-direction
};

class KineticArg : public Arg
{
  public:
    explicit KineticArg();
    explicit KineticArg(std::string);
    ~KineticArg();
    void Decode(std::string); // class function that imports parameters from a file
    DispersionType dispersion_type; // Kinetic operator dispersion relation
    Float mass; // mass for D operator
};


#endif

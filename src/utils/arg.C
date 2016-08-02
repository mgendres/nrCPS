#include <string>
#include <iostream>
using namespace std;
#include "arg.h"
#include "verbose.h"
#include "error.h"

Arg::Arg()
{
  const char* fname = "Arg::Arg()";
  VRB.Func(fname);
}

Arg::~Arg()
{
  const char* fname = "Arg::~Arg()";
  VRB.Func(fname);
}

void Arg::OpenFile(std::string filename)
{
  const char* fname = "void Arg::OpenFile(std::string)";
  VRB.Func(fname);
  file.open(filename.c_str());
  if (!file) { ERR.FileR(fname, filename.c_str()); }
}

void Arg::CloseFile()
{
  const char* fname = "void Arg::CloseFile()";
  VRB.Func(fname);
  file.close();
}

template <typename T>
int Arg::ReadLine(T &param) {
  const char* fname = "int Arg::ReadLine(T &)";
  VRB.Func(fname);
  file >> param;
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(bool &param)
{
  const char* fname = "int Arg::ReadLine(bool &)";
  VRB.Func(fname);
  file >> input;
  if (input=="1") {
    param=true;
  } else if (input=="0") {
    param=false;
  } else {
    ERR.General(fname, "Unable to recognize bool");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(RandomType &param)
{
  const char* fname = "int Arg::ReadLine(RandomType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="RANDOM_TYPE_RAN0") {
    param = RANDOM_TYPE_RAN0;
  } else if (input=="RANDOM_TYPE_RAN1") {
    param = RANDOM_TYPE_RAN1;
  } else if (input=="RANDOM_TYPE_RAN2") {
    param = RANDOM_TYPE_RAN2;
  } else if (input=="RANDOM_TYPE_RAN3") {
    param = RANDOM_TYPE_RAN3;
  } else if (input=="RANDOM_TYPE_RANLXS0") {
    param = RANDOM_TYPE_RANLXS0;
  } else if (input=="RANDOM_TYPE_RANLXS1") {
    param = RANDOM_TYPE_RANLXS1;
  } else if (input=="RANDOM_TYPE_RANLXS2") {
    param = RANDOM_TYPE_RANLXS2;
  } else if (input=="RANDOM_TYPE_RANLXD1") {
    param = RANDOM_TYPE_RANLXD1;
  } else if (input=="RANDOM_TYPE_RANLXD2") {
    param = RANDOM_TYPE_RANLXD2;
  } else {
    ERR.General(fname, "Unable to recognize RandomType");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(SeedType &param)
{
  const char* fname = "int Arg::ReadLine(SeedType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="SEED_TYPE_DEFAULT") {
    param = SEED_TYPE_DEFAULT;
  } else if (input=="SEED_TYPE_INPUT") {
    param = SEED_TYPE_INPUT;
  } else if (input=="SEED_TYPE_TIME") {
    param = SEED_TYPE_TIME;
  } else if (input=="SEED_TYPE_FILE") {
    param = SEED_TYPE_FILE;
  } else {
    ERR.General(fname, "Unable to recognize SeedType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(FieldType &param)
{
  const char* fname = "int Arg::ReadLine(FieldType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="FIELD_TYPE_GAUSS") {
    param = FIELD_TYPE_GAUSS;
  } else if (input=="FIELD_TYPE_COMPLEXGAUSS") {
    param = FIELD_TYPE_COMPLEXGAUSS;
  } else if (input=="FIELD_TYPE_ZTWO") {
    param = FIELD_TYPE_ZTWO;
  } else if (input=="FIELD_TYPE_ZTHREE") {
    param = FIELD_TYPE_ZTHREE;
  } else {
    ERR.General(fname, "Unable to recognize FieldType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(SourceType &param)
{
  const char* fname = "int Arg::ReadLine(SourceType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="SOURCE_TYPE_MOM") {
    param = SOURCE_TYPE_MOM;
  } else if (input=="SOURCE_TYPE_SHO") {
    param = SOURCE_TYPE_SHO;
  } else {
    ERR.General(fname, "Unable to recognize SourceType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(WavefuncType &param)
{
  const char* fname = "int Arg::ReadLine(WavefuncType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="WAVEFUNC_TYPE_NONE") {
    param = WAVEFUNC_TYPE_NONE;
  } else if (input=="WAVEFUNC_TYPE_UNIFORM") {
    param = WAVEFUNC_TYPE_UNIFORM;
  } else if (input=="WAVEFUNC_TYPE_GND") {
    param = WAVEFUNC_TYPE_GND;
  } else if (input=="WAVEFUNC_TYPE_TGND") {
    param = WAVEFUNC_TYPE_TGND;
  } else if (input=="WAVEFUNC_TYPE_PAIR1") {
    param = WAVEFUNC_TYPE_PAIR1;
  } else if (input=="WAVEFUNC_TYPE_PAIR2") {
    param = WAVEFUNC_TYPE_PAIR2;
  } else {
    ERR.General(fname, "Unable to recognize WavefuncType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(BoundaryType &param)
{
  const char* fname = "int Arg::ReadLine(BoundaryType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="BOUNDARY_TYPE_PRD") {
    param = BOUNDARY_TYPE_PRD;
  } else if (input=="BOUNDARY_TYPE_APRD") {
    param = BOUNDARY_TYPE_APRD;
  } else {
    ERR.General(fname, "Unable to recognize BoundaryType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(DispersionType &param)
{
  const char* fname = "int Arg::ReadLine(DispersionType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="DISPERSION_TYPE_STANDARD") {
    param = DISPERSION_TYPE_STANDARD;
  } else if (input=="DISPERSION_TYPE_QUADRATIC") {
    param = DISPERSION_TYPE_QUADRATIC;
  } else if (input=="DISPERSION_TYPE_PERFECT") {
    param = DISPERSION_TYPE_PERFECT;
  } else {
    ERR.General(fname, "Unable to recognize DispersionType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(PotentialForm &param)
{
  const char* fname = "int Arg::ReadLine(PotentialForm &)";
  VRB.Func(fname);
  file >> input;
  if (input=="POTENTIAL_FORM_NONE") {
    param = POTENTIAL_FORM_NONE;
  } else if (input=="POTENTIAL_FORM_HARMONIC") {
    param = POTENTIAL_FORM_HARMONIC;
  } else if (input=="POTENTIAL_FORM_COULOMB") {
    param = POTENTIAL_FORM_COULOMB;
  } else {
    ERR.General(fname, "Unable to recognize PotentialForm.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(PotentialType &param)
{
  const char* fname = "int Arg::ReadLine(PotentialType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="POTENTIAL_TYPE_LIN") {
    param = POTENTIAL_TYPE_LIN;
  } else if (input=="POTENTIAL_TYPE_EXP") {
    param = POTENTIAL_TYPE_EXP;
  } else {
    ERR.General(fname, "Unable to recognize PotentialType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(InteractionType &param)
{
  const char* fname = "int Arg::ReadLine(InteractionType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="INTERACTION_TYPE_NONE") {
    param = INTERACTION_TYPE_NONE;
  } else if  (input=="INTERACTION_TYPE_ONEMINUSXIINVSQ") {
    param = INTERACTION_TYPE_ONEMINUSXIINVSQ;
  } else if (input=="INTERACTION_TYPE_XISQMINUSONE") {
    param = INTERACTION_TYPE_XISQMINUSONE;
  } else {
    ERR.General(fname, "Unable to recognize InteractionType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(CutoffType &param)
{
  const char* fname = "int Arg::ReadLine(CutoffType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="CUTOFF_TYPE_NONE") {
    param = CUTOFF_TYPE_NONE;
  } else if  (input=="CUTOFF_TYPE_HARD") {
    param = CUTOFF_TYPE_HARD;
  } else {
    ERR.General(fname, "Unable to recognize CutoffType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(IoType &param)
{
  const char* fname = "int Arg::ReadLine(IoType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="IO_TYPE_ROOT") {
    param = IO_TYPE_ROOT;
  } else if  (input=="IO_TYPE_ALL") {
    param = IO_TYPE_ALL;
  } else {
    ERR.General(fname, "Unable to recognize IoType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

int Arg::ReadLine(DetType &param)
{
  const char* fname = "int Arg::ReadLine(DetType &)";
  VRB.Func(fname);
  file >> input;
  if (input=="DET_TYPE_STANDARD") {
    param = DET_TYPE_STANDARD;
  } else if  (input=="DET_TYPE_LOG") {
    param = DET_TYPE_LOG;
  } else {
    ERR.General(fname, "Unable to recognize DetType.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}


int Arg::ReadLine(Vector3 &param)
{
  const char* fname = "int Arg::ReadLine(Vector3 &)";
  VRB.Func(fname);

  int input;
  file >> input;
  param.SetX(input);
  file >> input;
  param.SetY(input);
  file >> input;
  param.SetZ(input);
  if ( file.eof() ) { return 1; } else { return 0; };

}

int Arg::ReadLine(OctahedralRep &param)
{
  const char* fname = "int Arg::ReadLine(OctahedralRep &)";
  VRB.Func(fname);
  file >> input;
  if        (input=="OCTAHEDRAL_REP_A1g") {
    param = OCTAHEDRAL_REP_A1g;
  } else if (input=="OCTAHEDRAL_REP_A2g") {
    param = OCTAHEDRAL_REP_A2g;
  } else if (input=="OCTAHEDRAL_REP_Eg") {
    param = OCTAHEDRAL_REP_Eg;
  } else if (input=="OCTAHEDRAL_REP_T1g") {
    param = OCTAHEDRAL_REP_T1g;
  } else if (input=="OCTAHEDRAL_REP_T2g") {
    param = OCTAHEDRAL_REP_T2g;
  } else if (input=="OCTAHEDRAL_REP_A1u") {
    param = OCTAHEDRAL_REP_A1u;
  } else if (input=="OCTAHEDRAL_REP_A2u") {
    param = OCTAHEDRAL_REP_A2u;
  } else if (input=="OCTAHEDRAL_REP_Eu") {
    param = OCTAHEDRAL_REP_Eu;
  } else if (input=="OCTAHEDRAL_REP_T1u") {
    param = OCTAHEDRAL_REP_T1u;
  } else if (input=="OCTAHEDRAL_REP_T2u") {
    param = OCTAHEDRAL_REP_T2u;
  } else {
    ERR.General(fname, "Unable to recognize OctahedralRep.");
  }
  if ( file.eof() ) { return 1; } else { return 0; };
}

//// Derived Arg classes 

DoArg::DoArg()
: t_sites(4),
  x_sites(4),
  y_sites(4),
  z_sites(4),
  x_boundary_type(BOUNDARY_TYPE_PRD),
  y_boundary_type(BOUNDARY_TYPE_PRD),
  z_boundary_type(BOUNDARY_TYPE_PRD),
  cutoff(1.0)
{
  const char* fname = "DoArg::DoArg()";
  VRB.Func(fname);
}

DoArg::DoArg(std::string filename)
{
  const char* fname = "DoArg::DoArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void DoArg::Decode(std::string filename)
{
  const char* fname = "void DoArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(t_sites);
  ReadLine(x_sites);
  ReadLine(y_sites);
  ReadLine(z_sites);
  ReadLine(x_boundary_type);
  ReadLine(y_boundary_type);
  ReadLine(z_boundary_type);
  ReadLine(cutoff);
  CloseFile();
}

DoArg::~DoArg()
{
  const char* fname = "DoArg::~DoArg()";
  VRB.Func(fname);
}

LatticeArg::LatticeArg()
: field_type(FIELD_TYPE_GAUSS)
{
  const char* fname = "LatticeArg::LatticeArg()";
  VRB.Func(fname);
}

LatticeArg::LatticeArg(std::string filename)
{
  const char* fname = "LatticeArg::LatticArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void LatticeArg::Decode(std::string filename)
{

  const char* fname = "void LatticeArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(field_type);
  CloseFile();
}

LatticeArg::~LatticeArg()
{
  const char* fname = "LatticeArg::~LatticeArg()";
  VRB.Func(fname);
}


RandomArg::RandomArg()
: random_type(RANDOM_TYPE_RAN2),
  seed_type(SEED_TYPE_INPUT),
  seed(7012),
  file_stem("rng")
{
  const char* fname = "RandomArg::RandomArg()";
  VRB.Func(fname);
}

RandomArg::RandomArg(std::string filename)
{
  const char* fname = "RandomArg::RandomArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void RandomArg::Decode(std::string filename)
{
  const char* fname = "void RandomArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(random_type);
  ReadLine(seed_type);
  ReadLine(seed);
  ReadLine(file_stem);
  CloseFile();
}

RandomArg::~RandomArg()
{
  const char* fname = "RandomArg::~RandomArg()";
  VRB.Func(fname);
}

PropagatorArg::PropagatorArg()
: n_fermions(1),
  file_stem("propagator")
{
  const char* fname = "PropagatorArg::PropagatorArg()";
  VRB.Func(fname);
}

PropagatorArg::PropagatorArg(std::string filename)
{
  const char* fname = "PropagatorArg::PropagatorArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void PropagatorArg::Decode(std::string filename)
{

  const char* fname = "void PropagatorArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(n_fermions);
  ReadLine(file_stem);
  CloseFile();
}

PropagatorArg::~PropagatorArg()
{
  const char* fname = "PropagatorArg::~PropagatorArg()";
  VRB.Func(fname);
}



OneBodyArg::OneBodyArg()
: source_type(SOURCE_TYPE_MOM),
  lambda1(1.0),
  lambda2(1.0),
  lambda3(1.0)
{
  const char* fname = "OneBodyArg::OneBodyArg()";
  VRB.Func(fname);
}

void OneBodyArg::Decode(std::string filename)
{

  const char* fname = "void OneBodyArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(source_type);
  ReadLine(lambda1);
  ReadLine(lambda2);
  ReadLine(lambda3);
  CloseFile();
}

OneBodyArg::OneBodyArg(std::string filename)
{
  const char* fname = "OneBodyArg::OneBodyArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

OneBodyArg::~OneBodyArg()
{
  const char* fname = "OneBodyArg::~OneBodyArg()";
  VRB.Func(fname);
}

TwoBodyArg::TwoBodyArg()
: wavefunc_type(WAVEFUNC_TYPE_UNIFORM),
  dispersion_type1(DISPERSION_TYPE_STANDARD),
  dispersion_type2(DISPERSION_TYPE_STANDARD),
  mass1(1.0),
  mass2(1.0),
  lambda(1.0),
  file_stem("two_body")
{
  const char* fname = "TwoBodyArg::TwoBodyArg()";
  VRB.Func(fname);
}

TwoBodyArg::TwoBodyArg(std::string filename)
{
  const char* fname = "TwoBodyArg::TwoBodyArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void TwoBodyArg::Decode(std::string filename)
{

  const char* fname = "void TwoArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(wavefunc_type);
  ReadLine(dispersion_type1);
  ReadLine(dispersion_type2);
  ReadLine(mass1);
  ReadLine(mass2);
  ReadLine(lambda);
  ReadLine(file_stem);
  CloseFile();
}

TwoBodyArg::~TwoBodyArg()
{
  const char* fname = "TwoBodyArg::~TwoBodyArg()";
  VRB.Func(fname);
}

GeneralizedTwoBodyArg::GeneralizedTwoBodyArg()
: wavefunc_type(WAVEFUNC_TYPE_UNIFORM),
//  rep(OCTAHEDRAL_REP_A1g),
  distinct_vec(0,0,0),
  lambda(1.0),
  file_stem("two_body")
{
  const char* fname = "GeneralizedTwoBodyArg::GeneralizedTwoBodyArg()";
  VRB.Func(fname);
}

GeneralizedTwoBodyArg::GeneralizedTwoBodyArg(std::string filename)
{
  const char* fname = "GeneralizedTwoBodyArg::GeneralizedTwoBodyArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void GeneralizedTwoBodyArg::Decode(std::string filename)
{

  const char* fname = "void GeneralizedTwoBodyArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(wavefunc_type);
//  ReadLine(rep);
  ReadLine(distinct_vec);
  ReadLine(lambda);
  ReadLine(file_stem);
  CloseFile();
}

GeneralizedTwoBodyArg::~GeneralizedTwoBodyArg()
{
  const char* fname = "GeneralizedTwoBodyArg::~GeneralizedTwoBodyArg()";
  VRB.Func(fname);
}

VerboseArg::VerboseArg()
: func_level(true),
  warn_level(true),
  result_level(true),
  flow_level(true),
  debug_level(true)
{
  const char* fname = "VerboseArg::VerboseArg()";
  VRB.Func(fname);
}

VerboseArg::VerboseArg(std::string filename)
{
  const char* fname = "VerbodyArg::VerboseArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void VerboseArg::Decode(std::string filename)
{

  const char* fname = "void VerboseArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(func_level);
  ReadLine(warn_level);
  ReadLine(result_level);
  ReadLine(flow_level);
  ReadLine(debug_level);
  CloseFile();
}

VerboseArg::~VerboseArg()
{
  const char* fname = "VerboseArg::~VerboseArg()";
  VRB.Func(fname);
}

EvoArg::EvoArg()
: start(0),
  unload_period(1),
  configurations(1)
{
  const char* fname = "EvoArg::EvoArg()";
  VRB.Func(fname);
}

EvoArg::EvoArg(std::string filename)
{
  const char* fname = "EvoArg::EvoArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void EvoArg::Decode(std::string filename)
{

  const char* fname = "void EvoArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(start);
  ReadLine(unload_period);
  ReadLine(configurations);
  CloseFile();
}

EvoArg::~EvoArg()
{
  const char* fname = "EvoArg::~EvoArg()";
  VRB.Func(fname);
}

MomentaArg::MomentaArg()
: dispersion_type(DISPERSION_TYPE_PERFECT),
  fermi_energy(1.0),
  mass(1.0)
{
  const char* fname = "MomentaArg::MomentaArg()";
  VRB.Func(fname);
}

MomentaArg::MomentaArg(std::string filename)
{
  const char* fname = "MomentaArg::MomentaArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void MomentaArg::Decode(std::string filename)
{

  const char* fname = "void MomentaArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(dispersion_type);
  ReadLine(fermi_energy);
  ReadLine(mass);
  CloseFile();
}

MomentaArg::~MomentaArg()
{
  const char* fname = "MomentaArg::~MomentaArg()";
  VRB.Func(fname);
}

InteractionArg::InteractionArg()
: mass(1.0),
  interaction_type(INTERACTION_TYPE_NONE),
  num_couplings(0)
{
  const char* fname = "InteractionArg::InteractionArg()";
  VRB.Func(fname);
}

InteractionArg::InteractionArg(std::string filename)
{
  const char* fname = "InteractionArg::InteractionArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void InteractionArg::Decode(std::string filename)
{

  const char* fname = "void InteractionArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(mass);
  ReadLine(interaction_type);
  Float coupling;
  while (true) {
    if ( ReadLine(coupling) ) break;
    couplings.push_back(coupling);
  }
  CloseFile();

  num_couplings = couplings.size();
  for (int n=0; n<num_couplings; ++n) {
    VRB.Flow(fname,"Coupling %d of %d = %e", n+1, num_couplings, couplings[n]);
  }

}

InteractionArg::~InteractionArg()
{
  const char* fname = "InteractionArg::~InteractionArg()";
  VRB.Func(fname);
}

PotentialArg::PotentialArg()
: potential_form(POTENTIAL_FORM_NONE),
  potential_type(POTENTIAL_TYPE_EXP),
  spring_constant1(0.0),
  spring_constant2(0.0),
  spring_constant3(0.0)
{
  const char* fname = "PotentialArg::PotentialArg()";
  VRB.Func(fname);
}

PotentialArg::PotentialArg(std::string filename)
{
  const char* fname = "PotentialArg::PotentialArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void PotentialArg::Decode(std::string filename)
{

  const char* fname = "void PotentialArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(potential_form);
  ReadLine(potential_type);
  ReadLine(spring_constant1);
  ReadLine(spring_constant2);
  ReadLine(spring_constant3);
  CloseFile();
}

PotentialArg::~PotentialArg()
{
  const char* fname = "PotentialArg::~PotentialArg()";
  VRB.Func(fname);
}

KineticArg::KineticArg()
: dispersion_type(DISPERSION_TYPE_STANDARD),
  mass(1.0)
{
  const char* fname = "KineticArg::KineticArg()";
  VRB.Func(fname);
}

KineticArg::KineticArg(std::string filename)
{
  const char* fname = "KineticArg::KineticArg(std::string)";
  VRB.Func(fname);
  Decode(filename);
}

void KineticArg::Decode(std::string filename)
{

  const char* fname = "void KineticArg::Decode(std::string)";
  VRB.Func(fname);

  OpenFile(filename);
  ReadLine(dispersion_type);
  ReadLine(mass);
  CloseFile();

}

KineticArg::~KineticArg()
{
  const char* fname = "KineticArg::~KineticArg()";
  VRB.Func(fname);
}


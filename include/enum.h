#include "config.h"

#ifndef INCLUDED_ENUM
#define INCLUDED_ENUM

enum RandomType
{
  RANDOM_TYPE_RAN0,
  RANDOM_TYPE_RAN1,
  RANDOM_TYPE_RAN2,
  RANDOM_TYPE_RAN3,
  RANDOM_TYPE_RANLXS0,
  RANDOM_TYPE_RANLXS1,
  RANDOM_TYPE_RANLXS2,
  RANDOM_TYPE_RANLXD1,
  RANDOM_TYPE_RANLXD2
};

enum SeedType
{
  SEED_TYPE_DEFAULT,
  SEED_TYPE_INPUT,
  SEED_TYPE_TIME,
  SEED_TYPE_FILE
}; 

enum FieldType
{
  FIELD_TYPE_GAUSS,
  FIELD_TYPE_COMPLEXGAUSS,
  FIELD_TYPE_ZTWO,
  FIELD_TYPE_ZTHREE
};

enum SourceType
{
  SOURCE_TYPE_MOM,
  SOURCE_TYPE_SHO
};

enum WavefuncType
{
  WAVEFUNC_TYPE_NONE,
  WAVEFUNC_TYPE_UNIFORM,
  WAVEFUNC_TYPE_GND,
  WAVEFUNC_TYPE_TGND,
  WAVEFUNC_TYPE_PAIR1,
  WAVEFUNC_TYPE_PAIR2
};

enum BoundaryType
{
  BOUNDARY_TYPE_PRD=0, // These must not be changed
  BOUNDARY_TYPE_APRD=1
};

enum DispersionType
{
  DISPERSION_TYPE_STANDARD,
  DISPERSION_TYPE_QUADRATIC,
  DISPERSION_TYPE_PERFECT
};

enum CutoffType
{
  CUTOFF_TYPE_NONE,
  CUTOFF_TYPE_HARD
};

enum PotentialType
{
  POTENTIAL_TYPE_LIN,
  POTENTIAL_TYPE_EXP
};

enum PotentialForm
{
  POTENTIAL_FORM_NONE,
  POTENTIAL_FORM_HARMONIC,
  POTENTIAL_FORM_COULOMB
};

enum DetType
{
  DET_TYPE_STANDARD,
  DET_TYPE_LOG
};

enum InteractionType
{
  INTERACTION_TYPE_NONE,
  INTERACTION_TYPE_XISQMINUSONE,
  INTERACTION_TYPE_ONEMINUSXIINVSQ
};

enum IoType
{
  IO_TYPE_ROOT,
  IO_TYPE_ALL
};

enum OctahedralRep
{
  OCTAHEDRAL_REP_A1g = 0, // These must not be changed
  OCTAHEDRAL_REP_A2g = 1,
  OCTAHEDRAL_REP_Eg  = 2,
  OCTAHEDRAL_REP_T1g = 3,
  OCTAHEDRAL_REP_T2g = 4,
  OCTAHEDRAL_REP_A1u = 5,
  OCTAHEDRAL_REP_A2u = 6,
  OCTAHEDRAL_REP_Eu  = 7,
  OCTAHEDRAL_REP_T1u = 8,
  OCTAHEDRAL_REP_T2u = 9
};

enum SubshellTexture
{
  SUBSHELL_TEXTURE_000 = 0, // These must not be changed
  SUBSHELL_TEXTURE_00J = 1,
  SUBSHELL_TEXTURE_0JJ = 2,
  SUBSHELL_TEXTURE_JJJ = 3,
  SUBSHELL_TEXTURE_0JK = 4,
  SUBSHELL_TEXTURE_JJK = 5,
  SUBSHELL_TEXTURE_JKK = 6,
  SUBSHELL_TEXTURE_JKL = 7
};

enum Parity
{
  PARITY_POS = 0,
  PARITY_NEG = 1
};


#endif
#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>
using namespace std;
#include "constants.h"
#include "vector3.h"
#include "oh_group.h"
#include "error.h"
#include "eigensystem.h"

const int L = 12;
const Float M = 5.0;
const Float M0 = 5.0;
const Float CUTOFF = 0.99999*PI;
const Float CUTOFF_SQ = CUTOFF * CUTOFF;

const Float c[3] = { 0.0, 1.0, 0.0 }; const int Nop = 3;

const Vector3 dims(L);
const Vector3 offset(1-L/2);
const OctahedralRep rr = OCTAHEDRAL_REP_T1u;
const OctahedralRep r = OCTAHEDRAL_REP_T1u;
const Vector3 distinct_vecs[488] = {
                                     Vector3(0,0,0),
                                     Vector3(0,0,1),
                                     Vector3(0,1,1),
                                     Vector3(1,1,1),
                                     Vector3(0,0,2),
                                     Vector3(0,1,2),
                                     Vector3(1,1,2),
                                     Vector3(0,2,2),
                                     Vector3(0,0,3),
                                     Vector3(1,2,2),
                                     Vector3(0,1,3),
                                     Vector3(1,1,3),
                                     Vector3(2,2,2),
                                     Vector3(0,2,3),
                                     Vector3(1,2,3),
                                     Vector3(0,0,4),
                                     Vector3(0,1,4),
                                     Vector3(2,2,3),
                                     Vector3(1,1,4),
                                     Vector3(0,3,3),
                                     Vector3(1,3,3),
                                     Vector3(0,2,4),
                                     Vector3(1,2,4),
                                     Vector3(2,3,3),
                                     Vector3(2,2,4),
                                     Vector3(0,0,5),
                                     Vector3(0,3,4),
                                     Vector3(0,1,5),
                                     Vector3(1,3,4),
                                     Vector3(1,1,5),
                                     Vector3(3,3,3),
                                     Vector3(0,2,5),
                                     Vector3(2,3,4),
                                     Vector3(1,2,5),
                                     Vector3(0,4,4),
                                     Vector3(2,2,5),
                                     Vector3(1,4,4),
                                     Vector3(0,3,5),
                                     Vector3(3,3,4),
                                     Vector3(1,3,5),
                                     Vector3(0,0,6),
                                     Vector3(2,4,4),
                                     Vector3(0,1,6),
                                     Vector3(1,1,6),
                                     Vector3(2,3,5),
                                     Vector3(0,2,6),
                                     Vector3(1,2,6),
                                     Vector3(0,4,5),
                                     Vector3(3,4,4),
                                     Vector3(1,4,5),
                                     Vector3(3,3,5),
                                     Vector3(2,2,6),
                                     Vector3(0,3,6),
                                     Vector3(2,4,5),
                                     Vector3(1,3,6),
                                     Vector3(4,4,4),
                                     Vector3(0,0,7),
                                     Vector3(2,3,6),
                                     Vector3(0,1,7),
                                     Vector3(0,5,5),
                                     Vector3(3,4,5),
                                     Vector3(1,1,7),
                                     Vector3(1,5,5),
                                     Vector3(0,4,6),
                                     Vector3(0,2,7),
                                     Vector3(1,4,6),
                                     Vector3(1,2,7),
                                     Vector3(3,3,6),
                                     Vector3(2,5,5),
                                     Vector3(2,4,6),
                                     Vector3(2,2,7),
                                     Vector3(4,4,5),
                                     Vector3(0,3,7),
                                     Vector3(1,3,7),
                                     Vector3(3,5,5),
                                     Vector3(0,5,6),
                                     Vector3(3,4,6),
                                     Vector3(2,3,7),
                                     Vector3(1,5,6),
                                     Vector3(0,0,8),
                                     Vector3(0,1,8),
                                     Vector3(0,4,7),
                                     Vector3(2,5,6),
                                     Vector3(1,1,8),
                                     Vector3(1,4,7),
                                     Vector3(4,5,5),
                                     Vector3(3,3,7),
                                     Vector3(0,2,8),
                                     Vector3(4,4,6),
                                     Vector3(1,2,8),
                                     Vector3(2,4,7),
                                     Vector3(3,5,6),
                                     Vector3(2,2,8),
                                     Vector3(0,6,6),
                                     Vector3(0,3,8),
                                     Vector3(1,6,6),
                                     Vector3(1,3,8),
                                     Vector3(0,5,7),
                                     Vector3(3,4,7),
                                     Vector3(1,5,7),
                                     Vector3(5,5,5),
                                     Vector3(2,6,6),
                                     Vector3(2,3,8),
                                     Vector3(4,5,6),
                                     Vector3(2,5,7),
                                     Vector3(0,4,8),
                                     Vector3(0,0,9),
                                     Vector3(1,4,8),
                                     Vector3(4,4,7),
                                     Vector3(3,6,6),
                                     Vector3(0,1,9),
                                     Vector3(3,3,8),
                                     Vector3(1,1,9),
                                     Vector3(3,5,7),
                                     Vector3(2,4,8),
                                     Vector3(0,2,9),
                                     Vector3(0,6,7),
                                     Vector3(1,2,9),
                                     Vector3(1,6,7),
                                     Vector3(5,5,6),
                                     Vector3(4,6,6),
                                     Vector3(2,2,9),
                                     Vector3(0,5,8),
                                     Vector3(3,4,8),
                                     Vector3(2,6,7),
                                     Vector3(0,3,9),
                                     Vector3(1,5,8),
                                     Vector3(4,5,7),
                                     Vector3(1,3,9),
                                     Vector3(2,5,8),
                                     Vector3(2,3,9),
                                     Vector3(3,6,7),
                                     Vector3(4,4,8),
                                     Vector3(0,4,9),
                                     Vector3(5,6,6),
                                     Vector3(1,4,9),
                                     Vector3(3,5,8),
                                     Vector3(0,7,7),
                                     Vector3(3,3,9),
                                     Vector3(1,7,7),
                                     Vector3(5,5,7),
                                     Vector3(0,0,10),
                                     Vector3(0,6,8),
                                     Vector3(0,1,10),
                                     Vector3(2,4,9),
                                     Vector3(1,6,8),
                                     Vector3(4,6,7),
                                     Vector3(1,1,10),
                                     Vector3(2,7,7),
                                     Vector3(0,2,10),
                                     Vector3(2,6,8),
                                     Vector3(1,2,10),
                                     Vector3(4,5,8),
                                     Vector3(0,5,9),
                                     Vector3(3,4,9),
                                     Vector3(1,5,9),
                                     Vector3(3,7,7),
                                     Vector3(2,2,10),
                                     Vector3(6,6,6),
                                     Vector3(0,3,10),
                                     Vector3(3,6,8),
                                     Vector3(1,3,10),
                                     Vector3(2,5,9),
                                     Vector3(5,6,7),
                                     Vector3(2,3,10),
                                     Vector3(4,4,9),
                                     Vector3(0,7,8),
                                     Vector3(1,7,8),
                                     Vector3(5,5,8),
                                     Vector3(4,7,7),
                                     Vector3(3,5,9),
                                     Vector3(0,4,10),
                                     Vector3(4,6,8),
                                     Vector3(1,4,10),
                                     Vector3(0,6,9),
                                     Vector3(2,7,8),
                                     Vector3(3,3,10),
                                     Vector3(1,6,9),
                                     Vector3(2,4,10),
                                     Vector3(0,0,11),
                                     Vector3(2,6,9),
                                     Vector3(6,6,7),
                                     Vector3(0,1,11),
                                     Vector3(4,5,9),
                                     Vector3(3,7,8),
                                     Vector3(1,1,11),
                                     Vector3(5,7,7),
                                     Vector3(0,2,11),
                                     Vector3(0,5,10),
                                     Vector3(3,4,10),
                                     Vector3(5,6,8),
                                     Vector3(1,2,11),
                                     Vector3(1,5,10),
                                     Vector3(3,6,9),
                                     Vector3(0,8,8),
                                     Vector3(2,2,11),
                                     Vector3(2,5,10),
                                     Vector3(1,8,8),
                                     Vector3(4,7,8),
                                     Vector3(0,3,11),
                                     Vector3(0,7,9),
                                     Vector3(1,3,11),
                                     Vector3(1,7,9),
                                     Vector3(5,5,9),
                                     Vector3(4,4,10),
                                     Vector3(2,8,8),
                                     Vector3(4,6,9),
                                     Vector3(2,3,11),
                                     Vector3(3,5,10),
                                     Vector3(2,7,9),
                                     Vector3(6,7,7),
                                     Vector3(0,6,10),
                                     Vector3(6,6,8),
                                     Vector3(0,4,11),
                                     Vector3(1,6,10),
                                     Vector3(3,8,8),
                                     Vector3(1,4,11),
                                     Vector3(5,7,8),
                                     Vector3(3,3,11),
                                     Vector3(3,7,9),
                                     Vector3(2,6,10),
                                     Vector3(2,4,11),
                                     Vector3(4,5,10),
                                     Vector3(5,6,9),
                                     Vector3(0,0,12),
                                     Vector3(4,8,8),
                                     Vector3(0,1,12),
                                     Vector3(3,6,10),
                                     Vector3(0,8,9),
                                     Vector3(1,1,12),
                                     Vector3(0,5,11),
                                     Vector3(3,4,11),
                                     Vector3(1,8,9),
                                     Vector3(4,7,9),
                                     Vector3(1,5,11),
                                     Vector3(7,7,7),
                                     Vector3(0,2,12),
                                     Vector3(1,2,12),
                                     Vector3(0,7,10),
                                     Vector3(2,8,9),
                                     Vector3(6,7,8),
                                     Vector3(2,5,11),
                                     Vector3(1,7,10),
                                     Vector3(5,5,10),
                                     Vector3(2,2,12),
                                     Vector3(4,6,10),
                                     Vector3(0,3,12),
                                     Vector3(4,4,11),
                                     Vector3(2,7,10),
                                     Vector3(6,6,9),
                                     Vector3(5,8,8),
                                     Vector3(1,3,12),
                                     Vector3(3,8,9),
                                     Vector3(3,5,11),
                                     Vector3(5,7,9),
                                     Vector3(2,3,12),
                                     Vector3(0,6,11),
                                     Vector3(1,6,11),
                                     Vector3(3,7,10),
                                     Vector3(0,4,12),
                                     Vector3(1,4,12),
                                     Vector3(2,6,11),
                                     Vector3(5,6,10),
                                     Vector3(4,8,9),
                                     Vector3(3,3,12),
                                     Vector3(4,5,11),
                                     Vector3(0,9,9),
                                     Vector3(7,7,8),
                                     Vector3(1,9,9),
                                     Vector3(2,4,12),
                                     Vector3(0,8,10),
                                     Vector3(6,8,8),
                                     Vector3(1,8,10),
                                     Vector3(4,7,10),
                                     Vector3(3,6,11),
                                     Vector3(2,9,9),
                                     Vector3(6,7,9),
                                     Vector3(2,8,10),
                                     Vector3(0,0,13),
                                     Vector3(0,5,12),
                                     Vector3(3,4,12),
                                     Vector3(0,1,13),
                                     Vector3(1,5,12),
                                     Vector3(0,7,11),
                                     Vector3(5,8,9),
                                     Vector3(1,1,13),
                                     Vector3(1,7,11),
                                     Vector3(5,5,11),
                                     Vector3(3,9,9),
                                     Vector3(6,6,10),
                                     Vector3(0,2,13),
                                     Vector3(2,5,12),
                                     Vector3(4,6,11),
                                     Vector3(3,8,10),
                                     Vector3(1,2,13),
                                     Vector3(2,7,11),
                                     Vector3(5,7,10),
                                     Vector3(4,4,12),
                                     Vector3(2,2,13),
                                     Vector3(7,8,8),
                                     Vector3(0,3,13),
                                     Vector3(3,5,12),
                                     Vector3(4,9,9),
                                     Vector3(1,3,13),
                                     Vector3(3,7,11),
                                     Vector3(7,7,9),
                                     Vector3(0,6,12),
                                     Vector3(4,8,10),
                                     Vector3(1,6,12),
                                     Vector3(0,9,10),
                                     Vector3(6,8,9),
                                     Vector3(2,3,13),
                                     Vector3(5,6,11),
                                     Vector3(1,9,10),
                                     Vector3(2,6,12),
                                     Vector3(0,4,13),
                                     Vector3(4,5,12),
                                     Vector3(0,8,11),
                                     Vector3(2,9,10),
                                     Vector3(6,7,10),
                                     Vector3(1,4,13),
                                     Vector3(1,8,11),
                                     Vector3(4,7,11),
                                     Vector3(3,3,13),
                                     Vector3(5,9,9),
                                     Vector3(2,4,13),
                                     Vector3(3,6,12),
                                     Vector3(2,8,11),
                                     Vector3(5,8,10),
                                     Vector3(3,9,10),
                                     Vector3(8,8,8),
                                     Vector3(0,7,12),
                                     Vector3(6,6,11),
                                     Vector3(0,5,13),
                                     Vector3(3,4,13),
                                     Vector3(1,7,12),
                                     Vector3(5,5,12),
                                     Vector3(3,8,11),
                                     Vector3(7,8,9),
                                     Vector3(1,5,13),
                                     Vector3(5,7,11),
                                     Vector3(0,0,14),
                                     Vector3(4,6,12),
                                     Vector3(0,1,14),
                                     Vector3(2,7,12),
                                     Vector3(4,9,10),
                                     Vector3(1,1,14),
                                     Vector3(2,5,13),
                                     Vector3(7,7,10),
                                     Vector3(6,9,9),
                                     Vector3(0,2,14),
                                     Vector3(0,10,10),
                                     Vector3(6,8,10),
                                     Vector3(1,2,14),
                                     Vector3(4,4,13),
                                     Vector3(4,8,11),
                                     Vector3(1,10,10),
                                     Vector3(3,7,12),
                                     Vector3(0,9,11),
                                     Vector3(3,5,13),
                                     Vector3(1,9,11),
                                     Vector3(2,2,14),
                                     Vector3(2,10,10),
                                     Vector3(0,3,14),
                                     Vector3(0,6,13),
                                     Vector3(5,6,12),
                                     Vector3(1,3,14),
                                     Vector3(1,6,13),
                                     Vector3(2,9,11),
                                     Vector3(6,7,11),
                                     Vector3(5,9,10),
                                     Vector3(0,8,12),
                                     Vector3(2,3,14),
                                     Vector3(2,6,13),
                                     Vector3(1,8,12),
                                     Vector3(4,7,12),
                                     Vector3(3,10,10),
                                     Vector3(8,8,9),
                                     Vector3(4,5,13),
                                     Vector3(5,8,11),
                                     Vector3(3,9,11),
                                     Vector3(7,9,9),
                                     Vector3(0,4,14),
                                     Vector3(2,8,12),
                                     Vector3(1,4,14),
                                     Vector3(7,8,10),
                                     Vector3(3,3,14),
                                     Vector3(3,6,13),
                                     Vector3(2,4,14),
                                     Vector3(6,6,12),
                                     Vector3(4,10,10),
                                     Vector3(3,8,12),
                                     Vector3(6,9,10),
                                     Vector3(0,7,13),
                                     Vector3(5,7,12),
                                     Vector3(4,9,11),
                                     Vector3(1,7,13),
                                     Vector3(5,5,13),
                                     Vector3(7,7,11),
                                     Vector3(0,5,14),
                                     Vector3(3,4,14),
                                     Vector3(4,6,13),
                                     Vector3(0,10,11),
                                     Vector3(6,8,11),
                                     Vector3(1,5,14),
                                     Vector3(2,7,13),
                                     Vector3(1,10,11),
                                     Vector3(4,8,12),
                                     Vector3(0,0,15),
                                     Vector3(2,5,14),
                                     Vector3(0,9,12),
                                     Vector3(2,10,11),
                                     Vector3(5,10,10),
                                     Vector3(0,1,15),
                                     Vector3(1,9,12),
                                     Vector3(8,9,9),
                                     Vector3(1,1,15),
                                     Vector3(3,7,13),
                                     Vector3(5,9,11),
                                     Vector3(4,4,14),
                                     Vector3(8,8,10),
                                     Vector3(0,2,15),
                                     Vector3(2,9,12),
                                     Vector3(6,7,12),
                                     Vector3(1,2,15),
                                     Vector3(3,5,14),
                                     Vector3(5,6,13),
                                     Vector3(3,10,11),
                                     Vector3(7,9,10),
                                     Vector3(0,6,14),
                                     Vector3(2,2,15),
                                     Vector3(1,6,14),
                                     Vector3(0,8,13),
                                     Vector3(5,8,12),
                                     Vector3(0,3,15),
                                     Vector3(1,8,13),
                                     Vector3(4,7,13),
                                     Vector3(3,9,12),
                                     Vector3(7,8,11),
                                     Vector3(1,3,15),
                                     Vector3(2,6,14),
                                     Vector3(6,10,10),
                                     Vector3(4,5,14),
                                     Vector3(2,8,13),
                                     Vector3(4,10,11),
                                     Vector3(2,3,15),
                                     Vector3(6,9,11),
                                     Vector3(0,4,15),
                                     Vector3(3,6,14),
                                     Vector3(6,6,13),
                                     Vector3(4,9,12),
                                     Vector3(1,4,15),
                                     Vector3(3,8,13),
                                     Vector3(7,7,12),
                                     Vector3(0,11,11),
                                     Vector3(3,3,15),
                                     Vector3(5,7,13),
                                     Vector3(1,11,11),
                                     Vector3(9,9,9),
                                     Vector3(0,10,12),
                                     Vector3(6,8,12),
                                     Vector3(2,4,15),
                                     Vector3(0,7,14),
                                     Vector3(1,10,12),
                                     Vector3(8,9,10),
                                     Vector3(1,7,14),
                                     Vector3(5,5,14),
                                     Vector3(2,11,11),
                                     Vector3(5,10,11),
                                     Vector3(4,6,14),
                                     Vector3(2,10,12),
                                     Vector3(2,7,14),
                                     Vector3(4,8,13),
                                     Vector3(8,8,11),
                                     Vector3(7,10,10),
                                     Vector3(0,5,15),
                                     Vector3(3,4,15),
                                     Vector3(0,9,13),
                                     Vector3(5,9,12),
                                     Vector3(1,5,15),
                                     Vector3(1,9,13),
                                     Vector3(3,11,11),
                                     Vector3(7,9,11),
                                     Vector3(3,10,12),
                                     Vector3(2,5,15),
                                     Vector3(3,7,14),
                                     Vector3(2,9,13),
                                     Vector3(6,7,13)
                                  };

int GetVecListSize()
{
  switch (L) {
    case 2:  return 1;  break;
    case 4:  return 4;  break;
    case 6:  return 8;  break;
    case 8:  return 15; break;
    case 10: return 25; break;
    case 12: return 40; break;
//  case 14: return ; break;
    case 16: return 79; break;
//  case 18: return ; break;
//  case 20: return ; break;
//  case 22: return ; break;
//  case 24: return ; break;
//  case 26: return ; break;
//  case 28: return ; break;
//  case 30: return ; break;
    case 32: return 488; break;
    default: cout << "Choose a different value for L sucka!" << endl; exit(EXIT_FAILURE);
  }
}

inline Float Xi(Float psq)
{
  if ( psq <= CUTOFF_SQ ) { return exp( psq/(2.0*M) ); }
  return numeric_limits<Float>::infinity();
}

inline Float C(Float psq)
{

  if ( psq > CUTOFF_SQ ) { psq = CUTOFF_SQ; }
  Float PSQ = M0 * (1.0 - exp( -psq/M0) );

  Float ret_val = 0.0;
  Float Op=1.0;
  for (int n=0; n<Nop; ++n) {
    ret_val += c[n] * Op;
    Op *= PSQ;
  }
  ret_val *= (4.0*PI/M0);
  ret_val /= (Float)(L*L*L);

  return ret_val;

}

void PrintVector(Vector3& v) {
  cout << "(" << v.GetX() << " " << v.GetY() << " " << v.GetZ() << ")" << endl;
}

inline Vector3 Mod(const Vector3& n, const  Vector3& dims, const Vector3& offset)
{
  return ( n + 20*dims-offset ) % dims + offset;
}



inline Float T(Vector3& pp, Vector3& qq, Vector3& p, Vector3& q)
{

  Float ret_val = 0.0;

  const Float w = (TWOPI/L)*(TWOPI/L);

  Float pp_sq = pp.L2NormSquared() * w;
  Float qq_sq = qq.L2NormSquared() * w;

  Float p_sq  =  p.L2NormSquared() * w;
  Float q_sq  =  q.L2NormSquared() * w;

  Vector3 ppp = Mod(pp-p,dims,offset);
  Vector3 qqq = Mod(qq-q,dims,offset);
  Float ppp_sq  =  ppp.L2NormSquared() * w;
  Float qqq_sq  =  qqq.L2NormSquared() * w;

  ret_val += KroneckerDelta(ppp) * KroneckerDelta(qqq);
  ret_val += KroneckerDelta(Mod(ppp+qqq,dims,offset)) * C(qqq_sq);
  ret_val /= sqrt( Xi(pp_sq) * Xi(qq_sq) * Xi(p_sq) * Xi(q_sq) );

  return ret_val;

}


void PrintMatrix(Float* m, int len) {
  cout << endl;
  for (int i=0; i<len; ++i) {
    for (int j=0; j<len; ++j) {
      printf("%.15f ", m[i+j*len]);
    }
    cout << endl;
  }

}


int main()
{

  const char* fname = "int Main()";

  cout << "Here we go!\n\n";

  int vec_list_size = GetVecListSize();

  int tmat_size = 0;
  OhGroup** oh_group = new OhGroup*[vec_list_size];
  for (int b=0; b<vec_list_size; ++b) {
    oh_group[b] = new OhGroup( distinct_vecs[b] , 1);
    tmat_size += oh_group[b]->GetSubshellDimension();
  }

  // Allocate memory for the transfer matrix
  Float* tmat = new Float [tmat_size*tmat_size];
  Vector3 v_ii, v_jj, v_i, v_j;
  Vector3 Rv_ii, Rv_jj, Rv_i, Rv_j;

  // Now the fun part...
  cout << "Constructing the matrix..." << endl;
  Float* tmat_p = tmat;
  int AA, A;

  AA=0;
  for (int II=0; II<vec_list_size; ++II) {
    for (int ii=0; ii< oh_group[II]->GetSubshellDimension(); ++ii) {
      v_ii = oh_group[II]->GetVector(ii);      

      A=0;
      for (int  I=0;  I<vec_list_size;  ++I) {
        for (int  i=0;  i< oh_group[ I]->GetSubshellDimension();  ++i) {
          v_i  = oh_group[ I]->GetVector( i);      

          if (AA>=A) {
            OhGroup Rii(v_ii, 0); // In this case, can use 1 as well
            OhGroup Ri(v_i, 0);
            
            *tmat_p = 0.0;
            for (int mm=0; mm<Rii.GetSubshellDimension(); ++mm) {
              Rv_ii = Rii.GetVector(mm);
              Rv_jj = -1 * Rv_ii;
              for (int m=0; m<Ri.GetSubshellDimension(); ++m) {
                Rv_i = Ri.GetVector(m);
                Rv_j  = -1 * Rv_i;
                *tmat_p += Rii.GetWeight(rr,mm) * Ri.GetWeight(r,m) * T( Rv_ii, Rv_jj, Rv_i, Rv_j );
              }
            }

          }

          tmat_p++;
          A++;
        }
      }

      AA++;
    }
  }

//  PrintMatrix(tmat, tmat_size);

  if (1) {
    cout << "Diagonalizing..." << endl;
    // Okay... now compute some eigenvalues
    char jobz = 'N';
    char uplo = 'U';
    Float w[tmat_size];
    Eigensystem( &jobz, &uplo, &tmat_size, tmat, &tmat_size, w);

    for (int i = 0; i < tmat_size; ++i) {
      double eval_i = w[i];
      if (eval_i > 1e-8) printf ("eigenvalue = %.15g\n", eval_i);
    }
  }

  delete tmat;
  for (int b=0; b<vec_list_size; ++b) { delete oh_group[b]; }
  delete [] oh_group;

  return(EXIT_SUCCESS);

}



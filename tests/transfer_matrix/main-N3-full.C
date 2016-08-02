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

const Float M = 5.0;
const Float M0 = 5.0;
const Float CUTOFF = 0.99999*PI;
const Float CUTOFF_SQ = CUTOFF * CUTOFF;

//const int L = 4;
//const Float c[1] = { 0.6730676670288 }; const int Nop = 1;
//const Float c[1] = { 0.6730676670288 }; const int Nop = 1;
//const Float c[2] = { 0.3334771210495575, 0.15520549778861692 }; const int Nop = 2;

const int L = 6;
const Float c[1] = { 0.689183613936 }; const int Nop = 1;
//const Float c[2] = { 0.4280911579576602, 0.1128065028146968 }; const int Nop = 2;
//const Float c[3] = { 0.6713580119423893, -0.2382942408671079, 0.07753490483228229 }; const int Nop = 3;

const Vector3 dims(L);
const Vector3 offset(1-L/2);
const OctahedralRep rep = OCTAHEDRAL_REP_A1g;
const Vector3 distinct_vecs[40] = { Vector3(0, 0, 0),
                                    Vector3(0, 0, 1),
                                    Vector3(0, 1, 1),
                                    Vector3(1, 1, 1),
                                    Vector3(0, 0, 2),
                                    Vector3(0, 1, 2),
                                    Vector3(1, 1, 2),
                                    Vector3(0, 2, 2),
                                    Vector3(0, 0, 3),
                                    Vector3(1, 2, 2),
                                    Vector3(0, 1, 3),
                                    Vector3(1, 1, 3),
                                    Vector3(2, 2, 2),
                                    Vector3(0, 2, 3),
                                    Vector3(1, 2, 3),
                                    Vector3(0, 0, 4),
                                    Vector3(0, 1, 4),
                                    Vector3(2, 2, 3),
                                    Vector3(1, 1, 4),
                                    Vector3(0, 3, 3),
                                    Vector3(1, 3, 3),
                                    Vector3(0, 2, 4),
                                    Vector3(1, 2, 4),
                                    Vector3(2, 3, 3),
                                    Vector3(2, 2, 4),
                                    Vector3(0, 0, 5),
                                    Vector3(0, 3, 4),
                                    Vector3(0, 1, 5),
                                    Vector3(1, 3, 4),
                                    Vector3(1, 1, 5),
                                    Vector3(3, 3, 3),
                                    Vector3(0, 2, 5),
                                    Vector3(2, 3, 4),
                                    Vector3(1, 2, 5),
                                    Vector3(0, 4, 4),
                                    Vector3(2, 2, 5),
                                    Vector3(1, 4, 4),
                                    Vector3(0, 3, 5),
                                    Vector3(3, 3, 4),
                                    Vector3(1, 3, 5)
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
    default: cout << "Choose a different value for L sucka!" << endl; exit(EXIT_FAILURE);
  }
}




inline Float Xi(Float psq)
{
  if ( psq < CUTOFF_SQ ) { return exp( psq/(2.0*M) ); }
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

inline Vector3 Mod(const Vector3& n, const  Vector3& dims, const Vector3& offset)
{
  return ( n + 20*dims-offset ) % dims + offset;
}

inline Float T(Vector3& pp, Vector3& qq, Vector3& kk, Vector3& p, Vector3& q, Vector3& k)
{

  Float ret_val = 0.0;

  const Float w = (TWOPI/L)*(TWOPI/L);

  Float pp_sq = pp.L2NormSquared() * w;
  Float qq_sq = qq.L2NormSquared() * w;
  Float kk_sq = kk.L2NormSquared() * w;

  Float p_sq  =  p.L2NormSquared() * w;
  Float q_sq  =  q.L2NormSquared() * w;
  Float k_sq  =  k.L2NormSquared() * w;

  Vector3 ppp = Mod(pp-p,dims,offset);
  Vector3 qqq = Mod(qq-q,dims,offset);
  Vector3 kkk = Mod(kk-k,dims,offset);
  Float ppp_sq  =  ppp.L2NormSquared() * w;
  Float qqq_sq  =  qqq.L2NormSquared() * w;
  Float kkk_sq  =  kkk.L2NormSquared() * w;

  ret_val += KroneckerDelta(ppp) * KroneckerDelta(qqq) * KroneckerDelta(kkk);
  ret_val += KroneckerDelta(ppp) * KroneckerDelta( Mod( qqq+kkk, dims, offset) ) * C(kkk_sq);
  ret_val += KroneckerDelta(qqq) * KroneckerDelta( Mod( kkk+ppp, dims, offset) ) * C(ppp_sq);
  ret_val += KroneckerDelta(kkk) * KroneckerDelta( Mod( ppp+qqq, dims, offset) ) * C(qqq_sq);
  ret_val /= sqrt( Xi(pp_sq) * Xi(qq_sq) * Xi(kk_sq) * Xi(p_sq) * Xi(q_sq) * Xi(k_sq) );

  return ret_val;

}

inline Float T_uud_A(Vector3& pp, Vector3& qq, Vector3& kk, Vector3& p, Vector3& q, Vector3& k)
{
  Float ret_val=0.0;
  ret_val += T(pp,qq,kk,p,q,k);
  ret_val -= T(qq,pp,kk,p,q,k);
  ret_val -= T(pp,qq,kk,q,p,k);
  ret_val += T(qq,pp,kk,q,p,k);
  return ret_val/4.0;
}

void PrintVector(Vector3& v) {
  cout << "(" << v.GetX() << " " << v.GetY() << " " << v.GetZ() << ")" << endl;
}



void PrintMatrix(Float* m, int len) {
  cout << endl;
  for (int i=0; i<len; ++i) {
    for (int j=0; j<len; ++j) {
      printf("%2.6f ", m[i+j*len]);
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
  tmat_size *= tmat_size;
  cout << tmat_size << endl;

  // Allocate memory for the transfer matrix
  Float* tmat = new Float [tmat_size*tmat_size];
  Vector3 v_ii, v_jj, v_kk, v_i, v_j, v_k;

  // Now the fun part...
  cout << "Constructing the matrix..." << endl;
  Float* tmat_p = tmat;
  int AA, A;

  AA=0;
  for (int II=0; II<vec_list_size; ++II) 
  for (int JJ=0; JJ<vec_list_size; ++JJ) {

    for (int ii=0; ii< oh_group[II]->GetSubshellDimension(); ++ii) 
    for (int jj=0; jj< oh_group[JJ]->GetSubshellDimension(); ++jj) {
      v_ii = oh_group[II]->GetVector(ii);      
      v_jj = oh_group[JJ]->GetVector(jj);      
      v_kk = -1 * Mod( v_ii+v_jj ,dims, offset );

      A=0;
      for (int  I=0;  I<vec_list_size;  ++I) 
      for (int  J=0;  J<vec_list_size;  ++J) {

        for (int  i=0;  i< oh_group[ I]->GetSubshellDimension();  ++i) 
        for (int  j=0;  j< oh_group[ J]->GetSubshellDimension();  ++j) {
          v_i  = oh_group[ I]->GetVector( i);      
          v_j  = oh_group[ J]->GetVector( j);      
          v_k  = -1 * Mod( v_i+v_j , dims, offset );


          if (AA>=A) { *tmat_p =  T_uud_A( v_ii, v_jj, v_kk, v_i, v_j, v_k  ); }

          tmat_p++;

          A++;
        }
      }
      AA++;
    }
    //cout << II << " " << JJ << endl;
  }

  if (1) {
    cout << "Diagonalizing..." << endl;

    char jobz = 'N';
    char uplo = 'U';
    Float w[tmat_size];
    Eigensystem( &jobz, &uplo, &tmat_size, tmat, &tmat_size, w);

    for (int i = 0; i < tmat_size; ++i) {
      double eval_i = w[i];
      //if (eval_i > 1e-8) printf ("eigenvalue = %.15f\n", eval_i );
      if (eval_i > 1e-8) printf ("eigenvalue = %.15f\n", -log(eval_i) / ( (TWOPI/L)*(TWOPI/L) / (2.0*M)) );
    }
  }




  delete [] tmat;
  for (int b=0; b<vec_list_size; ++b) { delete oh_group[b]; }
  delete [] oh_group;

  return(EXIT_SUCCESS);

}

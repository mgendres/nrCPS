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

const int using_z2 = 1;
const int antisymmetrize = 1;
const int L = 4;
const Float M = 5.0;
const Float M0 = 5.0;
const Float CUTOFF = 0.99999*PI;
const Float CUTOFF_SQ = CUTOFF * CUTOFF;
const Float c[1] = { 0.6730676670288431 }; const int Nop = 1;
//const Float c[3] = { 0.33347712104954197, 0.15520549778862391, 0.0 }; const int Nop = 3;
//const Float c[3] = { 0.9656387836088537, -0.494142463854144, 0.08136513173330479 }; const int Nop = 3;
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


inline Float T(Vector3& pp, Vector3& qq, Vector3& kk, Vector3& ll, Vector3& p, Vector3& q, Vector3& k, Vector3& l)
{

  Float ret_val = 0.0;

  const Float w = (TWOPI/L)*(TWOPI/L);

  Float pp_sq = pp.L2NormSquared() * w;
  Float qq_sq = qq.L2NormSquared() * w;
  Float kk_sq = kk.L2NormSquared() * w;
  Float ll_sq = ll.L2NormSquared() * w;

  Float p_sq  =  p.L2NormSquared() * w;
  Float q_sq  =  q.L2NormSquared() * w;
  Float k_sq  =  k.L2NormSquared() * w;
  Float l_sq  =  l.L2NormSquared() * w;

  Vector3 ppp = Mod(pp-p,dims,offset);
  Vector3 qqq = Mod(qq-q,dims,offset);
  Vector3 kkk = Mod(kk-k,dims,offset);
  Vector3 lll = Mod(ll-l,dims,offset);

  Float ppp_sq  =  ppp.L2NormSquared() * w;
  Float qqq_sq  =  qqq.L2NormSquared() * w;
  Float kkk_sq  =  kkk.L2NormSquared() * w;
  Float lll_sq  =  lll.L2NormSquared() * w;

  ret_val += KroneckerDelta(ppp) * KroneckerDelta(qqq) * KroneckerDelta(kkk) * KroneckerDelta(ll,l);

  ret_val += KroneckerDelta(ppp) * KroneckerDelta(qqq) * KroneckerDelta( Mod(kkk+lll,dims,offset) ) * C(kkk_sq);
  ret_val += KroneckerDelta(ppp) * KroneckerDelta(kkk) * KroneckerDelta( Mod(qqq+lll,dims,offset) ) * C(lll_sq);
  ret_val += KroneckerDelta(ppp) * KroneckerDelta(lll) * KroneckerDelta( Mod(kkk+qqq,dims,offset) ) * C(kkk_sq);
  ret_val += KroneckerDelta(qqq) * KroneckerDelta(kkk) * KroneckerDelta( Mod(ppp+lll,dims,offset) ) * C(ppp_sq);
  ret_val += KroneckerDelta(qqq) * KroneckerDelta(lll) * KroneckerDelta( Mod(kkk+ppp,dims,offset) ) * C(kkk_sq);
  ret_val += KroneckerDelta(kkk) * KroneckerDelta(lll) * KroneckerDelta( Mod(qqq+ppp,dims,offset) ) * C(qqq_sq);

  ret_val += KroneckerDelta( Mod(ppp+qqq,dims,offset) ) * KroneckerDelta( Mod(kkk+lll,dims,offset) ) * C(ppp_sq) * C(kkk_sq);
  ret_val += KroneckerDelta( Mod(ppp+kkk,dims,offset) ) * KroneckerDelta( Mod(qqq+lll,dims,offset) ) * C(ppp_sq) * C(lll_sq);
  ret_val += KroneckerDelta( Mod(ppp+lll,dims,offset) ) * KroneckerDelta( Mod(qqq+kkk,dims,offset) ) * C(ppp_sq) * C(qqq_sq);

  // This term vanishes for fermions, but is present for bosons when Z2 fields are used
  // This stems from the fact that < \phi^4 > = 1 (for Z2) rather than 3 (for Gaussian)
  // Note that this term is 1/V suppressed since there is one fewer delta functions.
  if (using_z2) {
    ret_val -= 2.0 * KroneckerDelta( Mod(ppp+qqq+kkk+lll, dims, offset) ) \
               * sqrt( C(ppp_sq) * C(qqq_sq) * C(kkk_sq) * C(lll_sq) ) / (L*L*L);
  }

  ret_val /= sqrt( Xi(pp_sq) * Xi(qq_sq) * Xi(kk_sq) * Xi(ll_sq) * Xi(p_sq) * Xi(q_sq) * Xi(k_sq) * Xi(l_sq) );

  return ret_val;

}


inline Float T_uudd_A(Vector3& pp, Vector3& qq, Vector3& kk, Vector3& ll, Vector3& p, Vector3& q, Vector3& k, Vector3& l)
{
  Float ret_val=0.0;

  ret_val += T(pp,qq,kk,ll,p,q,k,l);
  ret_val -= T(qq,pp,kk,ll,p,q,k,l);
  ret_val -= T(pp,qq,ll,kk,p,q,k,l);
  ret_val += T(qq,pp,ll,kk,p,q,k,l);

  ret_val -= T(pp,qq,kk,ll,q,p,k,l);
  ret_val += T(qq,pp,kk,ll,q,p,k,l);
  ret_val += T(pp,qq,ll,kk,q,p,k,l);
  ret_val -= T(qq,pp,ll,kk,q,p,k,l);

  ret_val -= T(pp,qq,kk,ll,p,q,l,k);
  ret_val += T(qq,pp,kk,ll,p,q,l,k);
  ret_val += T(pp,qq,ll,kk,p,q,l,k);
  ret_val -= T(qq,pp,ll,kk,p,q,l,k);

  ret_val += T(pp,qq,kk,ll,q,p,l,k);
  ret_val -= T(qq,pp,kk,ll,q,p,l,k);
  ret_val -= T(pp,qq,ll,kk,q,p,l,k);
  ret_val += T(qq,pp,ll,kk,q,p,l,k);

  return ret_val/16.0;

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

  int size = 0;
  OhGroup** oh_group = new OhGroup*[vec_list_size];
  for (int b=0; b<vec_list_size; ++b) {
    oh_group[b] = new OhGroup( distinct_vecs[b] , 1);
    size += oh_group[b]->GetSubshellDimension();
  }
  int tmat_size = size*size*size;
  cout << tmat_size << endl;

  // Allocate memory for the transfer matrix
  Float* tmat = new Float [tmat_size*tmat_size];
  Vector3 v_ii, v_jj, v_kk, v_ll, v_i, v_j, v_k, v_l;

  // Now the fun part...
  cout << "Constructing the matrix..." << endl;
  Float* tmat_p = tmat;
  int AA, A;

  AA=0;
  for (int II=0; II<vec_list_size; ++II)
  for (int JJ=0; JJ<vec_list_size; ++JJ)
  for (int KK=0; KK<vec_list_size; ++KK) {


    for (int ii=0; ii<oh_group[II]->GetSubshellDimension(); ++ii)
    for (int jj=0; jj<oh_group[JJ]->GetSubshellDimension(); ++jj)
    for (int kk=0; kk<oh_group[KK]->GetSubshellDimension(); ++kk) {
      v_ii = oh_group[II]->GetVector(ii);      
      v_jj = oh_group[JJ]->GetVector(jj);      
      v_kk = oh_group[KK]->GetVector(kk);      
      v_ll = -1 * Mod( v_ii+v_jj+v_kk , dims, offset);

      A=0;
      for (int  I=0;  I<vec_list_size;  ++I)
      for (int  J=0;  J<vec_list_size;  ++J)
      for (int  K=0;  K<vec_list_size;  ++K) {

        for (int  i=0;  i<oh_group[ I]->GetSubshellDimension();  ++i)
        for (int  j=0;  j<oh_group[ J]->GetSubshellDimension();  ++j)
        for (int  k=0;  k<oh_group[ K]->GetSubshellDimension();  ++k) {
          v_i  = oh_group[ I]->GetVector( i);      
          v_j  = oh_group[ J]->GetVector( j);      
          v_k  = oh_group[ K]->GetVector( k);      
          v_l  = -1 * Mod( v_i+v_j+v_k , dims, offset);

          if (AA>=A) {
            if (antisymmetrize) { 
              *tmat_p =  T_uudd_A( v_ii, v_jj, v_kk, v_ll, v_i, v_j, v_k, v_l  );
            } else {
              *tmat_p =  T( v_ii, v_jj, v_kk, v_ll, v_i, v_j, v_k, v_l  );
            }
          }

          tmat_p++;

          A++;
        }
      }
      AA++;
    }
    //cout << II << " " << JJ << " " << KK << endl;
  }

  if (1) {
    cout << "Diagonalizing..." << endl;
    char jobz = 'N';
    char uplo = 'U';
    Float w[tmat_size];
    Eigensystem( &jobz, &uplo, &tmat_size, tmat, &tmat_size, w);
 
    for (int i = 0; i < tmat_size; ++i) {
      double eval_i = w[i];
      //if (eval_i > 1e-8) printf ("eigenvalue = %.15f\n", eval_i);
      if (eval_i > 1e-8) printf ("eigenvalue = %.15f\n", -log(eval_i) / ( (TWOPI/L)*(TWOPI/L) / (2.0*M)) );
    }
  }



  delete tmat;
  for (int b=0; b<vec_list_size; ++b) { delete oh_group[b]; }
  delete [] oh_group;

  return(EXIT_SUCCESS);

}

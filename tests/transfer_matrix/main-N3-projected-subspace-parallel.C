#include <stdlib.h>

#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
#include <limits>
using namespace std;
#include "verbose.h"
#include "sysfunc.h"
#include "sysio.h"
#include "constants.h"
#include "vector3.h"
#include "oh_group.h"
#include "error.h"
#include "eigensystem.h"

#ifdef USE_MPI
#include <mpi.h>
#endif



const Float CUTOFF = 0.99999*PI;
const Float CUTOFF_SQ = CUTOFF * CUTOFF;

////// M = 5, L scan for various Cs ////////

const Float M = 5.0;
const Float M0 = 5.0;

//const int L = 4;
//const Float c[1] = { 0.6730676670288 }; const int Nop = 1;
//const Float c[2] = { 0.3334771210495575, 0.15520549778861692 }; const int Nop = 2;

//const int L = 6;
//const Float c[1] = { 0.689183613936 }; const int Nop = 1;
//const Float c[2] = { 0.4280911579576602, 0.1128065028146968 }; const int Nop = 2;
//const Float c[3] = { 0.6713580119423893, -0.2382942408671079, 0.07753490483228229 }; const int Nop = 3;


//const int L = 8;
//const Float c[1] = { 0.6809710120087866 }; const int Nop = 1;
//const Float c[2] = { 0.44287321879668906, 0.09807270183445489 }; const int Nop = 2;
//const Float c[3] = { 0.5394309658575215, -0.07573147731546012, 0.04139667178276734 }; const int Nop = 3;
//const Float c[4] = { 0.4872592878690856, 0.29804274431806, -0.21167497672063026, 0.04053113284110629 }; const int Nop = 4;

//const int L = 10;
//const Float c[1] = { 0.6848579405685634 }; const int Nop = 1;
//const Float c[2] = { 0.4552892362410875, 0.09394240839070421 }; const int Nop = 2;
//const Float c[3] = { 0.5176720644112826, -0.035281255606636676, 0.03203241265607874 }; const int Nop = 3;
//const Float c[4] = { 0.5268509968421444, 0.07601905719918768, -0.06692294420397041, 0.0178071796553536 }; const int Nop = 4;
//const Float c[5] = { 0.5852732481999289, -0.15077197954030375, 0.21209232853232207, -0.09741525985397437, 0.014552972725461882 }; const int Nop = 5;

const int L = 12;
//const Float c[1] = { 0.6797868177917205 }; const int Nop = 1;
//const Float c[2] = { 0.4559761508343166, 0.09003897181223035 }; const int Nop = 2;
//const Float c[3] = { 0.4998410139818714, -0.016061669654245422, 0.027058912890651685 }; const int Nop = 3;
const Float c[4] = { 0.514076019175471, 0.03555051390459175, -0.031524405983668444, 0.01143432801498804 }; const int Nop = 4;
//const Float c[5] = { 0.5440641387397142, -0.035487557967145776, 0.07704545156657798, -0.04331934366611181, 0.007838854957073964 }; const int Nop = 5;



////// 1C, M scan  ////////

//const Float M = 6.5;
//const Float M0 = 6.5;
//const int L = 8;
//const Float c[1] = { 0.6388736992559871 }; const int Nop = 1;
//const Float c[4] = { 0.4910342701956108, 0.15808037736052596, -0.10567388812869145, 0.018554955202901742 }; const int Nop = 4;

//const Float M = 10.0;
//const Float M0 = 10.0;
//const int L = 8;
//const Float c[1] = { 0.5909359520731624 }; const int Nop = 1;
//const Float c[4] = { 0.48555868880809927, 0.06592518210935944, -0.04493998578742857, 0.00720929440873537 }; const int Nop = 4;

//const Float M = 15.0;
//const Float M0 = 15.0;
//const int L = 8;
//const Float c[1] = { 0.5621577957667464 }; const int Nop = 1;
//const Float c[4] = { 0.48060172149992975, 0.029920643588353616, -0.02494430520962292, 0.003883920130807488 }; const int Nop = 4;

//const Float M = 20.0;
//const Float M0 = 20.0;
//const int L = 8;
//const Float c[1] = { 0.5480811357391961 }; const int Nop = 1;
//const Float c[4] = { 0.4781680544006946, 0.015696004739551146, -0.01794389779269297, 0.002806209004951569 }; const int Nop = 4;

//const Float M = 50.0;
//const Float M0 = 50.0;
//const int L = 8;
//const Float c[1] = { 0.5233321018117741 }; const int Nop = 1;
//const Float c[4] = { 0.4741307677079939, -0.005284049680845345, -0.008806127149472176, 0.0015003329005557938 }; const int Nop = 4;



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


void PrintMatrix(Float* m, int rlen, int clen) {
  for (int i=0; i<rlen; ++i) {
    for (int j=0; j<clen; ++j) {
      printf("%2.6f ", m[j+i*clen]);
    }
    cout << endl;
  }
  cout << endl;
  cout << endl;

}



void PrintMatrix(Float* m, int len) {
  cout << endl;
  cout << endl;
  for (int i=0; i<len; ++i) {
    for (int j=0; j<len; ++j) {
      printf("%2.6f ", m[j+i*len]);
    }
    cout << endl;
  }

}


int main()
{

  const char* fname = "int main()";
  const char* prog_name = "two_body";

  VRB.Flow(fname, "Running production code %s, compiled on %s at %s. ", prog_name, __DATE__, __TIME__ );
  VRB.LibInfo(fname);

  //---- Establish communications
  Comms::Initialize();
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  int vec_list_size = GetVecListSize();

  int tmat_size = 0;
  OhProjector** oh_projector = new OhProjector*[vec_list_size*vec_list_size];
  for (int I=0; I<vec_list_size; ++I) 
  for (int J=0; J<vec_list_size; ++J) {
    //cout << I << " " <<  J << endl;
//    if (Comms::Rank()==0) {
////      cout << Comms::Rank() << " " << Comms::Size() << endl;
//      cout << J+vec_list_size*I << " of " << vec_list_size*vec_list_size << endl;
//    }
    oh_projector[J+vec_list_size*I] = new OhProjector(distinct_vecs[I], distinct_vecs[J], rep);
    tmat_size += oh_projector[J+vec_list_size*I]->GetSubspaceDimension();
//    cout << tmat_size << endl;
  }
  cout << tmat_size << endl;

  OhGroup** oh_group = new OhGroup*[vec_list_size];
  for (int I=0; I<vec_list_size; ++I) {
    oh_group[I] = new OhGroup( distinct_vecs[I] , 1);
  }

  // Allocate memory for the transfer matrix
  Float* tmat = new Float [tmat_size*tmat_size];
  for (int i=0; i< tmat_size*tmat_size; ++i) { tmat[i]=0.0; }
  Vector3 v_ii, v_jj, v_kk, v_i, v_j, v_k;

  // Now the fun part...
  cout << "Constructing the matrix..." << endl;
  Float* tmat_p = tmat;

  int* block_dims = new int [vec_list_size*vec_list_size];

  Comms::Sync();
  for (int II=0; II<vec_list_size; ++II) 
  for (int JJ=0; JJ<vec_list_size; ++JJ) {

    //  cout << II << " " << JJ << endl;
    OhProjector* pp = oh_projector[JJ+vec_list_size*II];
    int PP_dim = pp->GetProjectorDimension();
    int PP_subdim = pp->GetSubspaceDimension();

    if (Comms::Rank()==0) {
//      cout << Comms::Rank() << " " << Comms::Size() << endl;
      cout << JJ+vec_list_size*II << " of " << vec_list_size*vec_list_size << endl;
    }

    block_dims[JJ+vec_list_size*II] = PP_subdim;
    if ( (JJ+vec_list_size*II)%Comms::Size() == Comms::Rank() )
    for (int  I=0;  I<vec_list_size;  ++I) 
    for (int  J=0;  J<vec_list_size;  ++J) {

        OhProjector* p = oh_projector[J+vec_list_size*I];
        int P_dim = p->GetProjectorDimension();
        int P_subdim = p->GetSubspaceDimension();

        // --allocate necessary memory for block
        Float* submat = new Float[PP_dim*P_dim];
        Float* submat_p = submat;

        // --fill in the matrix
        for (int ii=0; ii< oh_group[II]->GetSubshellDimension(); ++ii) 
        for (int jj=0; jj< oh_group[JJ]->GetSubshellDimension(); ++jj) {
          v_ii = oh_group[II]->GetVector(ii);      
          v_jj = oh_group[JJ]->GetVector(jj);      
          v_kk  = -1 * Mod( v_ii+v_jj , dims, offset );

          for (int  i=0;  i< oh_group[ I]->GetSubshellDimension();  ++i) 
          for (int  j=0;  j< oh_group[ J]->GetSubshellDimension();  ++j) {
            v_i  = oh_group[ I]->GetVector( i);      
            v_j  = oh_group[ J]->GetVector( j);      
            v_k  = -1 * Mod( v_i+v_j , dims, offset );

            *submat_p = T_uud_A( v_ii, v_jj, v_kk, v_i, v_j, v_k  );
            submat_p++;

          }
        }

        if (0) {
          cout << "**" <<  JJ+vec_list_size*II << " " << J+vec_list_size*I << endl;
          PrintMatrix(submat, PP_dim, P_dim);
        }

        // --perform projection
        for (int aa=0; aa<PP_subdim; ++aa) {
          Float* vex = pp->GetProjectorEigenvector(aa);
          for (int a=0; a<P_subdim; ++a) {
            Float* wex = p->GetProjectorEigenvector(a);

            Float inner = 0.0;
            for (int bb=0; bb<PP_dim; ++bb) {
              for (int b=0; b<P_dim; ++b) {
                inner += vex[bb] * submat[b+P_dim*bb] * wex[b];
              }
            }

            *tmat_p = inner;
            tmat_p++;
          }
          tmat_p += tmat_size - P_subdim;
        }
        tmat_p += P_subdim - PP_subdim*tmat_size;

        delete [] submat;

    }
    tmat_p += PP_subdim*tmat_size - tmat_size;
    if ( (JJ+vec_list_size*II)%Comms::Size() != Comms::Rank() ) {
      tmat_p += tmat_size;
    }
  }

//  for (int i=0; i< vec_list_size*vec_list_size; ++i) {
//    cout << block_dims[i] << endl;
//  }

//  for (int i=0; i<Comms::Size(); ++i) {
//    Comms::Sync();
//    if (Comms::Rank()==i ) {
//      PrintMatrix(tmat, tmat_size);
//      cout << endl;
//    }
//  }

  Comms::Sync();
  cout << "Moving data." << endl;
  MPI_Status status;
  tmat_p =tmat;
  for (int i=0; i<vec_list_size*vec_list_size; ++i) {
    Comms::Sync();
    if ( (i%Comms::Size()) == Comms::Rank() ) {
#ifdef USE_SINGLE
      MPI_Send( tmat_p, block_dims[i]*tmat_size, MPI_FLOAT, 0, 1,  MPI_COMM_WORLD );
#endif
#ifdef USE_DOUBLE
      MPI_Send( tmat_p, block_dims[i]*tmat_size, MPI_DOUBLE, 0, 1,  MPI_COMM_WORLD );
#endif
#ifdef USE_LONG_DOUBLE
      MPI_Send( tmat_p, block_dims[i]*tmat_size, MPI_LONG_DOUBLE, 0, 1,  MPI_COMM_WORLD );
#endif

    }
    if (Comms::Rank()==0) {
#ifdef USE_SINGLE
      MPI_Recv( tmat_p, block_dims[i]*tmat_size, MPI_FLOAT, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
#endif
#ifdef USE_DOUBLE
      MPI_Recv( tmat_p, block_dims[i]*tmat_size, MPI_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
#endif
#ifdef USE_LONG_DOUBLE
      MPI_Recv( tmat_p, block_dims[i]*tmat_size, MPI_LONG_DOUBLE, MPI_ANY_SOURCE, 1, MPI_COMM_WORLD, &status);
#endif

    }
    Comms::Sync();
    tmat_p += block_dims[i]*tmat_size;
  }

//  Comms::Sync();
//  if (Comms::Rank()==0) {
//    PrintMatrix(tmat, tmat_size);
//  }

  if (1&& (Comms::Rank()==0) ) {
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
  for (int I=0; I<vec_list_size; ++I) { delete oh_group[I]; }
  delete [] oh_group;
  for (int I=0; I<vec_list_size; ++I) 
  for (int J=0; J<vec_list_size; ++J) {
    delete oh_projector[J+vec_list_size*I];
  }
  delete [] oh_projector;
  delete [] block_dims;

  Comms::Finalize();
  return(EXIT_SUCCESS);

}

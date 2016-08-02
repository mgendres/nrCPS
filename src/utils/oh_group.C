#include <stdio.h>
#include <stdarg.h>
#include <iostream>
using namespace std;
#include "oh_group.h"
#include "global_job_parameter.h"
#include "error.h"
#include "eigensystem.h"

const int OhGroup::dim[10] = {1, 1, 2, 3, 3, 1, 1, 2, 3, 3};

const int OhGroup::chi[48][10] = { \
{  1,  1,  2,  3,  3,  1,  1,  2,  3,  3 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1, -1,  0,  0,  1,  1, -1,  0,  0 }, \
{  1,  1,  2, -1, -1,  1,  1,  2, -1, -1 }, \
{  1,  1,  2, -1, -1,  1,  1,  2, -1, -1 }, \
{  1,  1,  2, -1, -1,  1,  1,  2, -1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0,  1, -1,  1, -1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1,  1, -1,  0, -1,  1 }, \
{  1,  1,  2,  3,  3, -1, -1, -2, -3, -3 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1, -1,  0,  0, -1, -1,  1,  0,  0 }, \
{  1,  1,  2, -1, -1, -1, -1, -2,  1,  1 }, \
{  1,  1,  2, -1, -1, -1, -1, -2,  1,  1 }, \
{  1,  1,  2, -1, -1, -1, -1, -2,  1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0,  1, -1, -1,  1,  0, -1,  1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }, \
{  1, -1,  0, -1,  1, -1,  1,  0,  1, -1 }  \
};                                  

const int OhGroup::subshell_dimensions[8] = { 1, 6, 12, 8, 24, 24, 24, 48 };

OhGroup::OhGroup(const Vector3& vec, int reduced_flag)
: reduced(reduced_flag)
{

  const char* fname = "OhGroup::OhGroup(const Vector3&, int)";

  subshell_texture = GetSubshellTexture( GetDistinctVector(vec) );

  if (reduced) {
    subshell_dimension = subshell_dimensions[subshell_texture];
  } else {
    subshell_dimension = 48;
  }

  int x = vec.GetX();
  int y = vec.GetY();
  int z = vec.GetZ();

  // Generate all the possible vectors
  vectors[0]  =  Vector3( x,  y,  z); 
  vectors[1]  =  Vector3( z,  x,  y); 
  vectors[2]  =  Vector3( y, -z, -x); 
  vectors[3]  =  Vector3(-y, -z,  x); 
  vectors[4]  =  Vector3(-z, -x,  y); 
  vectors[5]  =  Vector3(-y,  z, -x); 
  vectors[6]  =  Vector3( z, -x, -y); 
  vectors[7]  =  Vector3(-z,  x, -y); 
  vectors[8]  =  Vector3( y,  z,  x); 
  vectors[9]  =  Vector3( x, -y, -z); 
  vectors[10] =  Vector3(-x,  y, -z); 
  vectors[11] =  Vector3(-x, -y,  z); 
  vectors[12] =  Vector3( x,  z, -y); 
  vectors[13] =  Vector3(-z,  y,  x); 
  vectors[14] =  Vector3( y, -x,  z); 
  vectors[15] =  Vector3( x, -z,  y); 
  vectors[16] =  Vector3( z,  y, -x); 
  vectors[17] =  Vector3(-y,  x,  z); 
  vectors[18] =  Vector3(-x,  z,  y); 
  vectors[19] =  Vector3( z, -y,  x); 
  vectors[20] =  Vector3( y,  x, -z); 
  vectors[21] =  Vector3(-x, -z, -y); 
  vectors[22] =  Vector3(-z, -y, -x); 
  vectors[23] =  Vector3(-y, -x, -z); 

  // And corresponding opposite parity 
  for (int k=0; k<24; ++k) { vectors[k+24] = -1*vectors[k]; }

  // Initialize the weights
  for (int j=0; j<48; ++j) {
    for (int k=0; k<10; ++k) {
      weights[j][k] = 0.0;
    }
  }

  if (reduced) {  // Computed wieghts for unique vectors

    // Identify vectors that are the same;
    // Assign an integer to each vector ranging from
    // 1 to subshell_dimension, stores in classifier array
    int classifier[48];
    for (int k=0; k<48; ++k) { classifier[k]=-1; }
  
    int count=0;
    for (int j=0; j<48; ++j) {
      if (classifier[j]==-1) {
        classifier[j]=count;
        for (int k=j; k<48; ++k) {
          if (classifier[k]==-1) {
            if (vectors[k]==vectors[j]) { classifier[k]=count; }
          }
        }
        count++;
      }
    }
  
    // Now that vectors are classified, lets sort them out and compute the weights

    // Begin search for vectors
    for (int j=0; j<subshell_dimension; ++j) {
      int setQ=0;
      for (int k=j; k<48; ++k) {
        if (classifier[k] == j) {
  
          // Set jth vector
          if (setQ==0) {
            vectors[j] = vectors[k];
            setQ=1;
            classifier[j] = j;
          }
  
          // sum up weights for each representation
          for (int l=0; l<10; ++l) {
            weights[j][l] += chi[k][l] * dim[l]/48.0;
          }
          
        }
      }
    }
  
 
  } else {  // Otherwise computed wieghts for all vectors

    for (int j=0; j<subshell_dimension; ++j) {
      for (int l=0; l<10; ++l) {
        weights[j][l] += chi[j][l] * dim[l]/48.0;
      }
    }

  }

  // Finally determine opposite parity indices
  opposite_parity_index = new int [subshell_dimension];
  for (int j=0; j<subshell_dimension; ++j) {
    for (int k=0; k<subshell_dimension; ++k) {
      if (vectors[j] == -1*vectors[k]) {
        opposite_parity_index[j] = k;
      }
    }
  }

}


OhGroup::~OhGroup()
{
  delete[] opposite_parity_index;
}

Vector3 OhGroup::GetVector(int k)
{
  const char* fname = "OhGroup::GetVector(int)";
  if ( (k<0) || (k>=subshell_dimension) ) { ERR.General(fname, "Argument out of bounds."); }
  return vectors[k];
}

Vector3 OhGroup::GetVector(int k, Parity parity)
{
  const char* fname = "OhGroup::GetVector(int, Parity)";
  if ( (k<0) || (k>=subshell_dimension) ) { ERR.General(fname, "Argument out of bounds."); }
  if (parity==PARITY_NEG) { return vectors[ opposite_parity_index[k] ]; }
  return vectors[k];
}

int OhGroup::GetOppositeParityIndex(int k)
{
  const char* fname = "OhGroup::GetOppositeParityIndex(int)";
  if ( (k<0) || (k>=subshell_dimension) ) { ERR.General(fname, "Argument out of bounds."); }
  return opposite_parity_index[k];
}

Float OhGroup::GetWeight(OctahedralRep rep, int k) {
  const char* fname = "OhGroup::GetWeight(OctahedralRep, int)";
  if ( (k<0) || (k>=subshell_dimension) ) { ERR.General(fname, "Argument out of bounds."); }
  return weights[k][rep];
}

int OhGroup::GetSubshellDimension() { return subshell_dimension; }

int OhGroup::GetSubshellTexture() { return subshell_texture; }


Vector3 OhGroup::GetDistinctVector(const Vector3& vec)
{
  int x = abs( vec.GetX() ); 
  int y = abs( vec.GetY() ); 
  int z = abs( vec.GetZ() ); 
  int w;
  if (x>y) { w=x; x=y; y=w; } 
  if (y>z) { w=y; y=z; z=w; } 
  if (x>y) { w=x; x=y; y=w; }
  return Vector3(x,y,z); 
}


SubshellTexture OhGroup::GetSubshellTexture(const Vector3& vec)
{

  Vector3 distinct_vec( GetDistinctVector(vec) );
  int x = distinct_vec.GetX(); 
  int y = distinct_vec.GetY(); 
  int z = distinct_vec.GetZ(); 

  if (x==0) {
    if (x==y) {
      if (y==z) {
        return SUBSHELL_TEXTURE_000;
      } else {
        return SUBSHELL_TEXTURE_00J;
      }
    } else {
      if (y==z) {
        return SUBSHELL_TEXTURE_0JJ;
      } else {
        return SUBSHELL_TEXTURE_0JK;
      }
    }
  } else {
    if (x==y) {
      if (y==z) {
        return SUBSHELL_TEXTURE_JJJ;
      } else {
        return SUBSHELL_TEXTURE_JJK;
      }
    } else {
      if (y==z) {
        return SUBSHELL_TEXTURE_JKK;
      } else {
        return SUBSHELL_TEXTURE_JKL;
      }
    }
  }

}


/////////////
/////////////


OhProjector::OhProjector(const Vector3& distinct_p, OctahedralRep rep)
{

  OhGroup p_vecs(distinct_p, 1);

  int p_dim = p_vecs.GetSubshellDimension();

  u_dim = p_dim;

  u = new Float[ u_dim * u_dim ];
  Vector3 p, pp, Rpp;
  Float elem;
  Float* u_p = u;

  for (int II=0; II < p_dim; ++II) {
    pp = p_vecs.GetVector(II);
    OhGroup Rpp_vecs(pp, 0); 

    for (int I=0; I < p_dim; ++I) {
      p = p_vecs.GetVector(I);

      elem = 0.0;
      for (int k=0; k<Rpp_vecs.GetSubshellDimension(); ++k) {
        Rpp = Rpp_vecs.GetVector(k);
        if ( p==Rpp ) { elem += Rpp_vecs.GetWeight(rep,k); }
      }

      *u_p = elem;
      u_p++;
    }
  }


  char jobz = 'V';
  char uplo = 'U';
  Float* w = new Float [u_dim];

  Eigensystem( &jobz, &uplo, &u_dim, u, &u_dim, w);

  subspace_dimension = 0;

  for (int i = 0; i < u_dim; ++i) {
    if ( (w[i]-1.0)*(w[i]-1.0) <=1e-6 ) { subspace_dimension++; }
  }

  delete [] w;


}


OhProjector::OhProjector(const Vector3& distinct_p, const Vector3& distinct_q, OctahedralRep rep)
{

  OhGroup p_vecs(distinct_p, 1);
  OhGroup q_vecs(distinct_q, 1);

  int p_dim = p_vecs.GetSubshellDimension();
  int q_dim = q_vecs.GetSubshellDimension();

  u_dim = p_dim*q_dim;

  u = new Float[ u_dim * u_dim ];
  Vector3 p, q, pp, qq, Rpp, Rqq;
  Float elem;
  Float* u_p = u;

  for (int II=0; II < p_dim; ++II) {
    pp = p_vecs.GetVector(II);
    OhGroup Rpp_vecs(pp, 0); 

    for (int JJ=0; JJ < q_dim; ++JJ) {
      qq = q_vecs.GetVector(JJ);
      OhGroup Rqq_vecs(qq, 0); 

      for (int I=0; I < p_dim; ++I) {
        p = p_vecs.GetVector(I);

        for (int J=0; J < q_dim; ++J) {
          q = q_vecs.GetVector(J);

          elem = 0.0;
          for (int k=0; k<Rpp_vecs.GetSubshellDimension(); ++k) {
            Rpp = Rpp_vecs.GetVector(k);
            Rqq = Rqq_vecs.GetVector(k);
            if ( (p==Rpp)&&(q==Rqq)) { elem += Rpp_vecs.GetWeight(rep,k); }
          }

          *u_p = elem;
          u_p++;
        }
      }
    }
  }

  char jobz = 'V';
  char uplo = 'U';
  Float* w = new Float [u_dim];

  Eigensystem( &jobz, &uplo, &u_dim, u, &u_dim, w);

  subspace_dimension = 0;
  for (int i = 0; i < u_dim; ++i) {
    if ( (w[i]-1.0)*(w[i]-1.0) <=1e-6 ) { subspace_dimension++; }
  }

  delete [] w;


}




OhProjector::OhProjector(const Vector3& distinct_p, const Vector3& distinct_q, const Vector3& distinct_r, OctahedralRep rep)
{

  OhGroup p_vecs(distinct_p, 1);
  OhGroup q_vecs(distinct_q, 1);
  OhGroup r_vecs(distinct_r, 1);

  int p_dim = p_vecs.GetSubshellDimension();
  int q_dim = q_vecs.GetSubshellDimension();
  int r_dim = r_vecs.GetSubshellDimension();

  u_dim = p_dim*q_dim*r_dim;

  u = new Float[ u_dim * u_dim ];
  Vector3 p, q, r, pp, qq, rr, Rpp, Rqq, Rrr;
  Float elem;
  Float* u_p = u;

  for (int II=0; II < p_dim; ++II) {
    pp = p_vecs.GetVector(II);
    OhGroup Rpp_vecs(pp, 0); 

    for (int JJ=0; JJ < q_dim; ++JJ) {
      qq = q_vecs.GetVector(JJ);
      OhGroup Rqq_vecs(qq, 0); 

      for (int KK=0; KK < r_dim; ++KK) {
        rr = r_vecs.GetVector(KK);
        OhGroup Rrr_vecs(rr, 0); 

        for (int I=0; I < p_dim; ++I) {
          p = p_vecs.GetVector(I);
        
          for (int J=0; J < q_dim; ++J) {
            q = q_vecs.GetVector(J);
        
            for (int K=0; K < r_dim; ++K) {
              r = r_vecs.GetVector(K);

              elem = 0.0;
              for (int g=0; g<Rpp_vecs.GetSubshellDimension(); ++g) {
                Rpp = Rpp_vecs.GetVector(g);
                Rqq = Rqq_vecs.GetVector(g);
                Rrr = Rrr_vecs.GetVector(g);
                if ( (p==Rpp)&&(q==Rqq)&&(r==Rrr)) { elem += Rpp_vecs.GetWeight(rep,g); }
              }
              
              *u_p = elem;
              u_p++;
            }
          }
        }
      }
    }
  }

  char jobz = 'V';
  char uplo = 'U';
  Float* w = new Float [u_dim];

  Eigensystem( &jobz, &uplo, &u_dim, u, &u_dim, w);

  subspace_dimension = 0;
  for (int i = 0; i < u_dim; ++i) {
    if ( (w[i]-1.0)*(w[i]-1.0) <=1e-6 ) { subspace_dimension++; }
  }

  delete [] w;


}

OhProjector::~OhProjector() { delete [] u; }
  
int OhProjector::GetProjectorDimension() { return u_dim; }

int OhProjector::GetSubspaceDimension() { return subspace_dimension; }

Float* OhProjector::GetProjectorEigenvector(int k)
{
  //const char* fname = "OhProjector::GetProjectorEigenvalue(int)";
  //if (k<0||k>=subspace_dimension) { ERR.General(fname, "Argument exceeds subspace dimension or is less than zero."); } 

  return u + u_dim*(u_dim-k-1);
}







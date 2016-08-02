#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "arg.h"
#include "global_job_parameter.h"
#include "vector3.h"
#include "oh_group.h"


int main()
{

  cout << "Here we go!\n\n";

  DoArg do_arg; // Container for lattice parameters
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg);

  Vector3 vecs[8];
  vecs[0].Set(0,0,0);
  vecs[1].Set(0,0,1);
  vecs[2].Set(0,1,1);
  vecs[3].Set(1,1,1);
  vecs[4].Set(0,1,2);
  vecs[5].Set(1,1,2);
  vecs[6].Set(1,2,2);
  vecs[7].Set(1,2,3);

  cout << endl;
  for (int i=0; i<8; ++i) {

    OhGroup oh_group( vecs[i] , 1);
    cout << "Distinct vector: ";
    cout << "(" << vecs[i].GetX() << " "<< vecs[i].GetY() << " "<< vecs[i].GetZ() << ")" << endl;
    cout << "Subshell dimension: " << oh_group.GetSubshellDimension() << ", Subshell texture: " << oh_group.GetSubshellTexture() << endl;
    Vector3 vec;
    for (int k=0; k<oh_group.GetSubshellDimension(); ++k) {
      vec = oh_group.GetVector(k) ;
      cout << "(" << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << ") : ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_A1g, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_A2g, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_Eg, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_T1g, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_T2g, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_A1u, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_A2u, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_Eu, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_T1u, k) << " ";
      cout << oh_group.GetWeight(OCTAHEDRAL_REP_T2u, k) << endl;
    }
    cout << endl;
    for (int k=0; k<oh_group.GetSubshellDimension(); ++k) {
      vec = oh_group.GetVector(k, PARITY_POS);
      cout << k << " : (" << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << ") : ";
      vec = oh_group.GetVector(k, PARITY_NEG);
      cout << oh_group.GetOppositeParityIndex(k) << " : (" << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << ") ";
      cout << endl;
    }
    cout << endl;
  }



  return(EXIT_SUCCESS);
}

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


int main()
{

  cout << "Here we go!\n\n";

  DoArg do_arg; // Container for lattice parameters
  do_arg.Decode("args/do.arg");
  GJP.Initialize(do_arg);

  cout << GJP.Xsites() << endl;
  cout << GJP.Ysites() << endl;
  cout << GJP.Zsites() << endl;

  Vector3 vec1(1,2,3);
  cout << "vect1: " << vec1.GetX() << " "<< vec1.GetY() << " "<< vec1.GetZ() << " " << vec1.L2NormSquared() << endl;

  Vector3 vec2(1,2,3);
  cout << "vec2: " << vec2.GetX() << " "<< vec2.GetY() << " "<< vec2.GetZ() << " " << vec2.L2NormSquared() << endl;

  vec2.SetX(7);
  vec2.SetY(8);
  vec2.SetZ(9);
  cout << "vec2: " << vec2.GetX() << " "<< vec2.GetY() << " "<< vec2.GetZ() << " " << vec2.L2NormSquared() << endl;

  vec2.Set(0, 1, 0);
  cout << "vec2: " << vec2.GetX() << " "<< vec2.GetY() << " "<< vec2.GetZ() << " " << vec2.L2NormSquared() << endl;

  Vector3 vec;
  cout << "vec: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared() << endl;
  vec = 1;
  cout << "vec: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared() << endl;
  vec = 0;
  cout << "vec: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared() << endl;

  vec = 0;
  vec += vec1;
  cout << "vec=0; vec+=vec1: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = vec2 + vec1;
  cout << "vec = vec1+vec2 : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = 0;
  vec -= vec1;
  cout << "vec=0; vec-=vec1: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = vec2 - vec1;
  cout << "vec = vec1-vec2 : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = vec1;
  vec *= 2;
  cout << "vec=vec1; vec*=2: " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = vec1 * 2;
  cout << "vec = vec1 * 2 : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  vec = 2 * vec1;
  cout << "vec = 2 * vec1 : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << " " << vec.L2NormSquared()<< endl;

  Vector3 dims = GJP.GetDimensions();
  for (int k=-20; k<20; ++k) {
    vec.Set(k,k,k);
    vec = IndexToVector3( Vector3ToIndex(vec, dims), dims );
    cout << k << ": " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << endl;
  }

  vec1.Set(1,2,3);
  vec2.Set(1,3,3);
  cout << (vec1==vec2) << endl;

  vec1=vec2;
  cout << (vec1==vec2) << endl;

  for (int k=-20; k<20; ++k) {
    vec.Set(k,k,k);
    vec = IndexToVector3( Vector3ToIndex(vec, dims), dims );
    cout << k << " : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ();
    vec = dims/2 - abs(vec - dims/2);
    cout << " : " << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << endl;
  }


  vec = Vector3(4,3,2)/Vector3(2,2,2);
  cout << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << endl;

  vec = 6/Vector3(2,2,2);
  cout << vec.GetX() << " "<< vec.GetY() << " "<< vec.GetZ() << endl;

  cout << KroneckerDelta(Vector3(1,3,3), Vector3(1,2,3)) << endl;
  cout << Dot(Vector3(1,3,3), Vector3(1,2,3)) << endl;

  return(EXIT_SUCCESS);
}

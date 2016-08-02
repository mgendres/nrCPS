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

const OctahedralRep r = OCTAHEDRAL_REP_T1u;
const Vector3 distinct_p(0,0,1);

int main()
{

  const char* fname = "int Main()";

  cout << "Here we go!\n\n";

  OhProjector oh_projector(distinct_p, r);
  cout << oh_projector.GetProjectorDimension() << endl;
  cout << oh_projector.GetSubspaceDimension() << endl;

  int projector_dimension = oh_projector.GetProjectorDimension();
  int subspace_dimension = oh_projector.GetSubspaceDimension();
  for (int i=0; i<subspace_dimension; ++i) {
    for (int j=0; j<subspace_dimension; ++j) {
      Float* eigv = oh_projector.GetProjectorEigenvector(i);
      Float* eigw = oh_projector.GetProjectorEigenvector(j);
      Float inner_product=0.0;
      for (int k=0; k<projector_dimension; ++k) {
        inner_product += eigv[k] * eigw[k];
      }
      printf("%.15f ", inner_product);
    }
    printf("\n");
  }

  return(EXIT_SUCCESS);

}

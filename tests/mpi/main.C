#include <stdlib.h>
#include <iostream>
#include <string>
#include <fstream>
#include <iomanip>
#include <complex>
using namespace std;
#include "sysfunc.h"
#include "sysio.h"

int main()
{

  const char* fname = "int main()";
  const char* prog_name = "two_body";

  //---- Establish communications
  Comms::Initialize();

  //---- Print something
  cout << "\nI'm " << Comms::Rank() << " of " << Comms::Size() << "\n";

  //---- Close communications
  Comms::Finalize();

  return(EXIT_SUCCESS);
}

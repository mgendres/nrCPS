#include <stdlib.h>
#include <iostream>
#include <ostream>
#include <limits>

using namespace std;

typedef std::numeric_limits< long double > ldl;
typedef std::numeric_limits< double > dl;
typedef std::numeric_limits< float > fl;

int main()
{

  cout << "long double:\n";
  cout << "\tdigits (bits):\t\t" << ldl::digits << endl;
  cout << "\tdigits (decimal):\t" << ldl::digits10 << endl;

  cout << endl;

  cout << "double:\n";
  cout << "\tdigits (bits):\t\t" << dl::digits << endl;
  cout << "\tdigits (decimal):\t" << dl::digits10 << endl;

  cout << endl;

  cout << "float:\n";
  cout << "\tdigits (bits):\t\t" << fl::digits << endl;
  cout << "\tdigits (decimal):\t" << fl::digits10 << endl;


  return(EXIT_SUCCESS);

}

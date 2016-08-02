#include <iostream>
#include <string>
using namespace std;
#include "cuda_utils.h"

#define N 10

template<typename T>
void printv(T array, int len)
{
  for (int i=0; i<len; ++i) { cout << array[i] << " "; }
  cout << endl;
}

int main()
{

  const char* fname = "int main()";
  cout << "Here we go!\n\n";

  Cuda::DeviceQ();

  return(EXIT_SUCCESS);

}

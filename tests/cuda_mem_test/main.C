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

  //---- Allocate arrays on host and device
  int a[N];
  int* dev_a;
  Cuda::Malloc( (void**)&dev_a, N*sizeof(int)   );

  //---- Initialize host array
  for (int i=0; i<N; ++i) { a[i] = i; }

  printv(a, N);

  //---- Load host array onto device 
  Cuda::MemCopy(dev_a, a, N*sizeof(int), cudaMemcpyHostToDevice); 

  //---- Clear host array
  for (int i=0; i<N; ++i) { a[i] = 0; }

  printv(a, N);

  //---- Load device array onto host
  Cuda::MemCopy(a, dev_a, N*sizeof(int), cudaMemcpyDeviceToHost);

  printv(a, N);

  //---- Free memory on device
  Cuda::Free( dev_a );

  return(EXIT_SUCCESS);

}

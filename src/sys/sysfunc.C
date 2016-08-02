#include <stdlib.h>
#include <sysfunc.h>
#include <cstddef>
#include <config.h>
#include "enum.h"
using namespace std;

#ifdef USE_MPI
#include <mpi.h>
#endif

//---- Some file scoped variables
namespace Comms
{
  static int rank = -1;
  static int size = -1;
}

bool Comms::initialized = false;

void Comms::Initialize()
{

  if (initialized) { return; }

#ifdef USE_MPI

  MPI_Init(NULL, NULL);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  MPI_Comm_size(MPI_COMM_WORLD, &size);

#else

  rank = 0;
  size = 1;

#endif

  initialized = true;

}

void Comms::Finalize() {

  if (!initialized) { return; }

#ifdef USE_MPI
  MPI_Finalize();
#else
#endif

}


void Comms::RaiseError() {
  Finalize();
  exit(EXIT_FAILURE);
}

int Comms::Rank()
{
  if (!initialized) { Initialize(); }
  return rank;
}

int Comms::Size()
{
  if (!initialized) { Initialize(); }
  return size;
}


void Comms::Sync() {

#ifdef USE_MPI
  if (!Comms::initialized) { Comms::Initialize(); }
  MPI_Barrier(MPI_COMM_WORLD);
#endif

}

complex<Float> Comms::GlobalSum(complex<Float> z)
{

  if (!initialized) { Initialize(); }

#ifdef USE_MPI
  Float sendbuf[2];
  Float recvbuf[2];

  sendbuf[0] = z.real();
  sendbuf[1] = z.imag();

  MPI_Barrier(MPI_COMM_WORLD);

#ifdef USE_SINGLE
  MPI_Reduce( sendbuf, recvbuf, 2, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );
#endif
#ifdef USE_DOUBLE
  MPI_Reduce( sendbuf, recvbuf, 2, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif
#ifdef USE_LONG_DOUBLE
  MPI_Reduce( sendbuf, recvbuf, 2, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif

  return complex<Float>(recvbuf[0], recvbuf[1]);
#else
  return z;
#endif

}


void Comms::GlobalSum(complex<Float>* z, int len)
{

  if (!initialized) { Initialize(); }

#ifdef USE_MPI
  Float sendbuf[2*len];
  Float recvbuf[2*len];

  int j, k;

  j = 0; k = 1;
  for (int i=0; i<len; ++i) {
    sendbuf[j] = z[i].real();
    sendbuf[k] = z[i].imag();
    j += 2; k += 2;
  }

  MPI_Barrier(MPI_COMM_WORLD);
#ifdef USE_SINGLE
  MPI_Reduce( sendbuf, recvbuf, 2*len, MPI_FLOAT, MPI_SUM, 0, MPI_COMM_WORLD );
#endif
#ifdef USE_DOUBLE
  MPI_Reduce( sendbuf, recvbuf, 2*len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif
#ifdef USE_LONG_DOUBLE
  MPI_Reduce( sendbuf, recvbuf, 2*len, MPI_LONG_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );
#endif

  j = 0; k = 1;
  for (int i=0; i<len; ++i) {
    z[i] = complex<Float>(recvbuf[j], recvbuf[k]);
    j += 2; k += 2;
  }

#endif

}

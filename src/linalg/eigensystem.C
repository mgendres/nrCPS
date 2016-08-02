#include <stdlib.h>
#include <stdio.h>
#include <Accelerate/Accelerate.h>
#include "eigensystem.h"

void Eigensystem( char* jobz, char* uplo, int* n, double* a, int* lda, double* w)
{

  /* Query and allocate the optimal workspace */
  int info=0;
  double wkopt;
  int lwork = -1;
  dsyev_( jobz, uplo, n, a, lda, w, &wkopt, &lwork, &info );

  lwork = (int)wkopt;
//  printf("%d\n", lwork);

  double* work = (double*)malloc( lwork*sizeof(double) );

  /* Solve eigenproblem */
  dsyev_( jobz, uplo, n, a, lda, w, work, &lwork, &info );

  /* Check for convergence */
  if( info > 0 ) { printf( "The algorithm failed to compute eigenvalues.\n" ); exit( 1 ); }

  /* Free workspace */
  free( (void*)work );


}

void Eigensystem( char* jobz, char* uplo, int* n, double* ap, double* w, int* ldz)
{

  if( *jobz =='V' ) { printf( "Sorry sucka! I'm too lazy to do this.\n" ); exit(1); }

  /* Query and allocate the optimal workspace */
  int info=0;

  double* work = (double*)malloc( 3 * (*n) * sizeof(double) );

  /* Solve eigenproblem */
  double* z; // Dummy
  dspev_( jobz, uplo, n, ap, w, z, ldz, work, &info );

  /* Check for convergence */
  if( info != 0 ) { printf( "The algorithm failed to compute eigenvalues.\n" ); exit( 1 ); }

  /* Free workspace */
  free( (void*)work );

}

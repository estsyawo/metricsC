/*
Invert matrices
*/

#include "solve.h"

void mat_inv( double *A, double *Ainv, int *n, char *solver )
{
int i,j;
double *X, *b;
// allocate memory for X and b
	X = malloc((*n)*sizeof(double));
	b = malloc((*n)*sizeof(double));

  for( j=0;j<*n;j++){ // run over columns of A
  // construct column of identity matrix
	for( i=0;i<*n;i++){
	if(i!=j){
    b[i] = 0.0 ;
    }else{
    b[i] = 1.0 ;
    } // end if
    X[i]=0.0; // initialise X to detect weird behaviour
	} // end for i construction b

  solve(A, b, X, n, solver); // solve for x in Ax=b

  for( i=0;i<*n;i++){ // extract x into Ainv
  Ainv[j*(*n) + i] = X[i] ;
  }//end for i extraction x

} // end for j

}


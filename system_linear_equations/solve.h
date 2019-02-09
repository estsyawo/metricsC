#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.

void gauss_seidel( double *A, double *b, double *phi, int *n);
void conjgrad(double *A, double *b, double *X, int *n);
void SOR( double *A, double *b, double *phi, int *n);

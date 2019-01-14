#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.
void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb);
void trans(double *a, double *atrans,int *nra,int *nca);
void gauss_seidel( double *A, double *b, double *phi, double *dev, int *n);
void conjgrad(double *A, double *b, double *X, int *n);
void SOR( double *A, double *b, double *phi, double *dev, int *n);

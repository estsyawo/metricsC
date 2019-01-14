/*header file for functions in matply.c*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.
void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb);
void trans(double *a, double *atrans,int *nra,int *nca);
double dotprod(double *a, double *b, int *n);
double max(double *x, int *n);
double min(double *x, int *n);

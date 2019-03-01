/*header file for functions in matply.c*/

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))
#define sgn(a) ((a) == 0 ? 0 : ((a) > 0  ? 1 : -1 )) // sign function
#define max_2(a,b) ((a) >= (b) ? (a) : (b))  // max of two values
#define min_2(a,b) ((a) <= (b) ? (a) : (b))  // min of two values


// declare function prototypes.
void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb);
void matply_sym(double *x, double *xx, int *nrx, int *ncx);
void matply_xty(double *x, double *y, double *xty, int *nrx, int *ncx, int *ncy);
void matply_xyt(double *x, double *y, double *xyt, int *nrx, int *ncx, int *nry);
void matply_sk1(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb, int *ica);
void trans(double *a, double *atrans,int *nra,int *nca);
double dotprod(double *a, double *b, int *n);
double dotprod_col_ex(double *a, double *b, int *nr, int *ica, int *icb);
void matply_ax(double *x, double *a, double *ax, int *n);
void vecadd(double *xa, double *xb, double *xm, int *n);
void vecsub(double *xa, double *xb, double *xm, int *n);
double max(double *x, int *n);
double min(double *x, int *n);
double norm_lp(double *x, int *n, int *p);
double norm_max(double *x, int *n);

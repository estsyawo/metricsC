#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#define absval(x) ((x) >=0.0 ? (x):(-(x)))

// declare function prototypes.
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_qrc(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_cord(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_qrMGS2(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx);
void linreg_qrMGS(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx);

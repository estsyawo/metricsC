/* 
Run linear regression for beta in Y = X*coefs + e. 
Solvers of system of normal equations available
 are 1. conjugate gradient, 2. gauss-seidel, and 3. SOR
*/

// prototype function declarations

#include "linreg.h"
#include "matply.h"

// linear regression with the conjugate gradient solver
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX,*tX,*XY;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = malloc(((*ncX)*(*ncX))*sizeof(double));
    XY = malloc((*nrX)*sizeof(double));
    tX = malloc(((*nrX)*(*ncX))*sizeof(double));
    
    trans(X, tX, nrX, ncX); // take X'
    matply(tX, X, XX, ncX, nrX, ncX); // take X'X
    matply(tX, Y, XY, ncX, nrX, &ncY); // take X'y
    
    conjgrad(XX, XY, coefs, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
    free(tX);
}

// linear regression with the gauss-seidel solver
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX, *tX, *XY, *dev;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = malloc(((*ncX)*(*ncX))*sizeof(double));
    XY = malloc((*nrX)*sizeof(double));
    tX = malloc(((*nrX)*(*ncX))*sizeof(double));
    dev = malloc((*nrX)*sizeof(double));
    
    trans(X, tX, nrX, ncX); // take X'
    matply(tX, X, XX, ncX, nrX, ncX); // take X'X
    matply(tX, Y, XY, ncX, nrX, &ncY); // take X'y
    
    gauss_seidel( XX, XY, coefs, dev, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
    free(tX);
}

// linear regression with the sor solver
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX, *tX, *XY, *dev;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = malloc(((*ncX)*(*ncX))*sizeof(double));
    XY = malloc((*nrX)*sizeof(double));
    tX = malloc(((*nrX)*(*ncX))*sizeof(double));
    dev = malloc((*nrX)*sizeof(double));
    
    trans(X, tX, nrX, ncX); // take X'
    matply(tX, X, XX, ncX, nrX, ncX); // take X'X
    matply(tX, Y, XY, ncX, nrX, &ncY); // take X'y
    
    SOR( XX, XY, coefs, dev, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
    free(tX);
}

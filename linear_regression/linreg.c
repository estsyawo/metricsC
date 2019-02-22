/* 
Run linear regression for beta in Y = X*coefs + e. 
Solvers of system of normal equations available
 are 1. conjugate gradient, 2. gauss-seidel, and 3. SOR
*/

#include "matply.h"
#include "solve.h"
#include "utils.h"

// linear regression with the conjugate gradient solver
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX,*tX,*XY;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = allocvector((*ncX)*(*ncX));
    XY = allocvector(*nrX);
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    
    conjgrad(XX, XY, coefs, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
}

// linear regression with the gauss-seidel solver
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX, *tX, *XY;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = allocvector((*ncX)*(*ncX));
    XY = allocvector(*nrX);

    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    
    gauss_seidel( XX, XY, coefs, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
}

// linear regression with the sor solver
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX, *tX, *XY, *dev;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = allocvector((*ncX)*(*ncX));
    XY = allocvector(*nrX);
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y
    
    SOR( XX, XY, coefs, ncX); // solve X'X*beta = X'Y
    
    // free allocated memory
    free(XX);
    free(XY);
}

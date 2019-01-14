//  linreg_test.c
//
//  Created by Emmanuel Tsyawo on 1/13/19.

/*
 Run linear regression for beta in Y = X*coefs + e.
 Solvers of system of normal equations available
 are 1. conjugate gradient, 2. gauss-seidel, and 3. SOR
 */
/*
 Compile using
 gcc -c linreg_test.c matply.c conjgrad.c gauss_seidel.c SOR.c
 gcc -o execLinreg_Test linreg_test.o matply.o conjgrad.o gauss_seidel.o SOR.o
 ./execLinreg_Test
 */


#include "linreg.h"
#include "matply.h"

// prototype function declarations
void linreg_cg(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_gs(double *Y, double *X, double *coefs, int *nrX, int *ncX);
void linreg_sor(double *Y, double *X, double *coefs, int *nrX, int *ncX);

// main function for code execution
int main( )
{
    int nrX = 10;
    int ncX = 2;
    double *coefs;
    
    coefs=malloc(ncX*sizeof(double));
    
    // design matrix X and response Y
    double X[20] = {1,1,1,1,1,1,1,1,1,1,
        -1.48,1.58,-0.96,-0.92,-2,-0.27,-0.32,-0.63,-0.11,0.43};
    double Y[10] = {2.48,-0.58,1.96,1.92,3,1.27,1.32,1.63,1.11,0.57};
    
    printf("The vector of true parameters is [1,-1].\n Calling linreg_cg( ) ... \n");
    linreg_cg(Y, X, coefs, &nrX, &ncX);
    printf("The solution is beta = [%.2f,%.2f]\n\n",coefs[0],coefs[1]);
    printf("The vector of true parameters is [1,-1].\n Calling linreg_gs( ) ... \n");
    linreg_gs(Y, X, coefs, &nrX, &ncX);
    printf("The solution is beta = [%.2f,%.2f]\n\n",coefs[0],coefs[1]);
    printf("The vector of true parameters is [1,-1].\n Calling linreg_sor( ) ... \n");
    linreg_sor(Y, X, coefs, &nrX, &ncX);
    printf("The solution is beta = [%.2f,%.2f]\n",coefs[0],coefs[1]);
    puts(" ");
}

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

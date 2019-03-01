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

// a wrapper for lsqqr - least squares via QR decomposition ||X*coefs-Y||^2
void linreg_qrc(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double *XX, *tX, *XY, *dev;
    int ncY;
    
    ncY = 1;
    // allocate memory
    XX = allocvector((*ncX)*(*ncX));
    XY = allocvector(*nrX);
    
    matply_sym(X, XX, nrX, ncX); // take X'X
    matply_xty(X, Y, XY, nrX, ncX, &ncY); // take X'Y

    lsqqr(XX, XY, coefs, ncX, ncX);
    
    // free allocated memory
    free(XX);
    free(XY);
}

// linear regression by coordinate descent

void linreg_cord(double *Y, double *X, double *coefs, int *nrX, int *ncX)
{
    double az = 1e-20, tol = 1e-10, *Xdot,*ervec,RSS0,RSS1,dev,dr;
    int ncY=1, k=0,j,id,maxk = 1000;
    // allocate memory
    Xdot = allocvector(*ncX); // to store sum of squares of columns in X
    ervec = allocvector(*nrX); // to store vector of length nrX
    // compute column-wise sum-of-squares
    //Xdot[0] = (double)(*nrX); // intercept's 1's
    for (j=1; j<=*ncX; j++) {
        Xdot[j-1]=dotprod_col_ex(X, X, nrX, &j, &j);
    }//end for
    
    // compute initial sum of squared of errors
    matply(X, coefs, ervec, nrX, ncX, &ncY); //compute ervec = X*coefs
    vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
    RSS0=dotprod(ervec, ervec, nrX); // compute initial residual sum of squares (RSS0)
    
    // commence iteration
    for(;;)
    {
        k+=1; //increment counter by 1
        for (j=0; j<*ncX; j++){
            // compute X*coefs excluding column j in X and element j in coefs
            matply_sk1(X, coefs,ervec, nrX, ncX, &ncY, &j); // ervec = X[,-j]*coef[-j]
            vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
            id = j+1; // dotprod_col_ex() needs id, that counts from 1,..,ncX
            dr=dotprod_col_ex(ervec, X, nrX, &ncY, &id); // compute dot(ervec,X[,j])
            coefs[j] = dr/Xdot[j]; // update coefs
        }// end for j
        
        // compute sum of squared residuals
        matply(X, coefs, ervec, nrX, ncX, &ncY); //compute ervec = X*coefs
        vecsub(Y, ervec, ervec, nrX); // compute ervec <-- Y-ervec
        RSS1=dotprod(ervec, ervec, nrX); // compute residual sum of squares (RSS1)
        
        dev = (RSS0-RSS1)/(az + RSS0); // compute deviation
        
        // check for convergence
        if (dev<=tol) {
            break;
        }
        
        if (k>=maxk) {
            printf("Warning: maximum number of iterations reached. \n");
            break;
        }
        RSS0 = RSS1; //update RSS0
    }// end for(;;)
    free(Xdot);
    free(ervec);
}


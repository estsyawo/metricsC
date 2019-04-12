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
        
        // compute deviation
        dev = (RSS0-RSS1)/(az + RSS0);
        
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
    // free allocated memory
    free(Xdot);
    free(ervec);
}

/***********************************************************************************/

// Reference: Algorithm 5.2.6 Matrix Computations, Golub & Van Loan 2013
// computes QR factorisation using the modified Gram-Schmidt algorithm
// A is overwritten by Q such that A = QR
// A is mxn, Q is mxn, R is nxn, m>=n
void qr_MGS(double *A, double *Q, double *R, int *m, int *n)
{
    int k,j,l,n1,n2,p=2;
    double sig;
    
    // begin main for loop
    for (k=0; k<*n; k++) {
        //n1 = k*(*m) + 1; n2 = k*(*m) + *m; // index k'th column in A
        sig = 0.0;
        for (l=0; l<*m; l++) {
            sig+= A[k*(*m)+l]*A[k*(*m)+l];
        }
        
        R[k*(*n)+k] = sqrt(sig); // compute R[k,k]
        for (l=0; l<*m; l++) {
            Q[k*(*m)+l] = A[k*(*m)+l]/R[k*(*n)+k]; //compute Q[,k]
        } // end for l
        if (k<(*n-1)) {
            for (j=(k+1); j<*n; j++) {
                sig=0.0;
                for (l=0; l<*m; l++) {
                    sig += Q[k*(*m)+l]*A[j*(*m)+l]; //compute Q
                } // end for l
                R[j*(*n)+k]=sig;
                for (l=0; l<*m; l++) {
                    A[j*(*m)+l] = A[j*(*m)+l] - Q[k*(*m)+l]*R[j*(*n)+k]; //compute Q
                } // end for l
                
            }// end for j
        } // end if
    }// end for k
    
} // end void function qr_MGS

/***********************************************************************************/
// this version overwrites A with Q
void qr_MGS2(double *A, double *R, int *m, int *n)
{
    int k,j,l,n1,n2,p=2;
    double *qv, sig;
    
    // allocate memory
    qv = allocvector(*m);
    
    // begin main for loop
    for (k=0; k<*n; k++) {
        sig = 0.0;
        for (l=0; l<*m; l++) {
            sig+= A[k*(*m)+l]*A[k*(*m)+l];
        }
        
        R[k*(*n)+k] = sqrt(sig); // compute R[k,k]
        
        for (l=0; l<*m; l++) {
            qv[l] = A[k*(*m)+l]/R[k*(*n)+k]; //compute Q[,k]
        } // end for l
        if (k<(*n-1)){
            for (j=(k+1); j<*n; j++) {
                sig=0.0;
                for (l=0; l<*m; l++) {
                    sig += qv[l]*A[j*(*m)+l]; //compute Q
                } // end for l
                R[j*(*n)+k]=sig;
                for (l=0; l<*m; l++) {
                    A[j*(*m)+l] = A[j*(*m)+l] - qv[l]*R[j*(*n)+k]; //compute Q
                } // end for l
                
            }// end for j
        } // end if
        // overwrite A[,k] with Q[,k]
        for (l=0; l<*m; l++) {
            A[k*(*m)+l] = qv[l]; //compute Q
        } // end for l
    }// end for k
    
} // end void function qr_MGS2
/***********************************************************************************/

// Reference: Algorithm 3.1.2 Matrix Computations, Golub & Van Loan 2013
// this function solves for x in Rx=b using a row-oriented back substitution. R is upper
// triangular. This algorithm overwrites b with the solution x.

void backsolve(double *R, double *b, int *n)
{
    int i,j;
    double sig;
    
    b[*n-1] = b[*n-1]/R[(*n-1)*(*n+1)];
    
    for (i=(*n-2); i>=0; i--) {
        sig = 0.0;
        for (j=(i+1); j<*n; j++) {
            sig+= (R[j*(*n)+i]*b[j]);
        }
        
        b[i] = (b[i] - sig)/R[i*(*n+1)];
    }
}

/***********************************************************************************/
//

void linreg_qrMGS2(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx)
{
    int ncy = 1;
    qr_MGS2(X, R, nrx, ncx); // QR decompose  X, X stores Q,
    matply_xty(X, Y, coefs, nrx, ncx, &ncy); // coefs <-- Q'Y
    backsolve(R, coefs, ncx); // coefs is overwritten with solution
}

void linreg_qrMGS(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx)
{
    int ncy = 1;
    double *Q;
    Q = allocvector((*nrx)*(*ncx)); // allocate memory
    qr_MGS(X, Q, R, nrx, ncx); // QR decompose  X, X stores Q,
    matply_xty(Q, Y, coefs, nrx, ncx, &ncy); // coefs <-- Q'Y
    backsolve(R, coefs, ncx); // coefs is overwritten with solution
}

/***********************************************************************************/
//  This file computes the QR factorisation of an mxn matrix A in QR where R is upper
//  triangular using the Modified Gram-Schmidt algorithm. It also computes the solves
//  the least squares problem.
//
//  Created by Emmanuel Tsyawo on 3/22/19.
//
#include "matply.h"
#include "utils.h"
/***********************************************************************************/
/* Compile using:
 gcc -c QR_MGS.c matply.c utils.c read_txt.c
 gcc -o exMGS QR_MGS.o matply.o utils.o read_txt.o
 ./exMGS
 */
/***********************************************************************************/
// prototype declaration
void qr_MGS(double *A, double *Q, double *R, int *m, int *n);
void qr_MGS2(double *A, double *R, int *m, int *n);
//void backsolve(double *R, double *b, int *n);
void linreg_qrMGS2(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx);
void linreg_qrMGS(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx);
// main programme
int main()
{
    int i,j, nr=12, nc = 11, nc1=nc-1;
    double *X, *X1, *Y, *dat, *R, *coefs;
    char *datname = "data.txt";
    
    Y = allocvector(nr); // allocate memory to Y
    X = allocvector(nr*nc1);
    X1 = allocvector(nr*nc1);
    R = allocvector(nc1*nc1);
    coefs = allocvector(nc1);
    // read in data
    dat=read_txt(datname, &nr, &nc );
    
    for (i=0; i<nr; i++) {
        for (j=0; j<(nc-1); j++) {
            X[j*nr+i]=dat[j*nr+i];
            Y[i]=dat[(nc-1)*nr+i];
        }
    }
    
    
    printf("Print out data set \n");
    printm(dat, &nr, &nc);
    
    printf("Print out matrix X \n");
    printm(X, &nr, &nc1);
    
    printf("Print out vector Y \n");
    printv(Y, &nr);
    
    qr_MGS(X, X1, R, &nr, &nc1);
    
    
    printf("Matrix Q is \n");
    printm(X1,&nr,&nc1);
    
    printf("The matrix R is \n");
    printm(R,&nc1,&nc1);
    
    printf("Can we recover matrix X?\n Compute QR =  ");
    matply(X1, R, X, &nr, &nc1, &nc1);
    printm(X,&nr,&nc1);
    
    printf("Evaluate Q'Q \n");
    matply_sym(X1, R, &nr, &nc1);
    printm(R,&nc1,&nc1);
    
    printf("Try version with overwriting \n");
    qr_MGS2(X, R, &nr, &nc1);
    
    printf("Matrix Q is \n");
    printm(X,&nr,&nc1);
    
    printf("The matrix R is \n");
    printm(R,&nc1,&nc1);
    
    
    printf("Solve the overdetermined system of equations problem ||X*coefs-Y||_2 \n");
    
    for (i=0; i<nr; i++) { //retrieve X again
        for (j=0; j<nc1; j++) {
            X[j*nr+i]=dat[j*nr+i];
        }
    }
    
    printf("First with version that does not overwrite X with Q\n");
    linreg_qrMGS(Y, X, coefs, R, &nr, &nc1);
    printf("The solution is \n");
    printv(coefs,&nc1);
    
    for (i=0; i<nr; i++) { //retrieve X again
        for (j=0; j<nc1; j++) {
            X[j*nr+i]=dat[j*nr+i];
        }
    }
    
    printf("Then with version that overwrites X with Q\n");
    linreg_qrMGS2(Y, X, coefs, R, &nr, &nc1);
    printf("The solution is \n");
    printv(coefs,&nc1);
    
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

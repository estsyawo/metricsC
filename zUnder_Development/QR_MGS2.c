/***********************************************************************************/
//  mat_decomp_MGS.c
//  
//
//  Created by Emmanuel Tsyawo on 3/22/19.
//
#include "matply.h"
#include "utils.h"
/***********************************************************************************/
/* Compile using:
 gcc -c QR_MGS2.c matply.c utils.c read_txt.c
 gcc -o exMGS2 QR_MGS2.o matply.o utils.o read_txt.o
 ./exMGS2
 */
/***********************************************************************************/
// prototype declaration
void qr_MGS(double *A, double *Q, double *R, int *m, int *n);
void qr_MGS2(double *A, double *R, int *m, int *n);
//void backsolve(double *R, double *b, int *n);
void linreg_qrMGS2(double *Y, double *X, double *coefs, double *R, int *nrx, int *ncx);

// main programme
int main()
{
    int nrX = 1000, ncX = 6, ny=1;
    double *coefs, *X, *Y, *R;
    char *datname;
    
    // allocate memory
    coefs=allocvector(ncX); Y = allocvector(nrX); R = allocvector(ncX*ncX);
    
    datname= "dat_lreg.txt"; // name of data set in folder
    // read in data
    X=read_txt(datname, &nrX, &ncX);
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nrX, &ncX);
    
    // true parameter values to be solved for
    double trcoefs[6]={1.2,0.5,1.0,0.0,1.5,-0.5};
    
    printf("True parameter values to be solved for");
    printm(trcoefs,&ny,&ncX);
    
    printf("Solve the overdetermined system of equations problem ||X*coefs-Y||_2 \n");
    
    linreg_qrMGS2(Y, X, coefs, R, &nrX, &ncX);
    printf("The solution is \n");
    printm(coefs,&ny,&ncX);

    
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

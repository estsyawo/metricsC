//
//  mat_decomp2.c
//  
//
//  Created by Emmanuel Tsyawo on 3/22/19.
/* Compile using:
 gcc -c mat_decomp2.c matply.c utils.c read_txt.c
 gcc -o exMD2 mat_decomp2.o matply.o utils.o read_txt.o
 ./exMD2
 */

#include "matply.h"
#include "utils.h"

// prototype declaration:
void house(double *A, double *v, double *beta, int *m, int *ia, int *ic);
void qr_household(double *A, int *m, int *n);
void back_accum(double *A, double *Q, int *m, int *n);
void linreg_hh(double *X, double *Y, int *nrx, int *ncx);

// main programme
int main()
{
    int i,j, nr=12, nc = 11, nc1=nc-1,itemp;
    double *X, *X1, *Y, *dat, *R, *coefs, *beta, *v,*Q;
    char *datname = "data.txt";
    
    Y = allocvector(nr); // allocate memory to Y
    X = allocvector(nr*nc1);
    Q = allocvector(nr*nr);
    v = allocvector(nr);
    R = allocvector(nc1*nc1);
    coefs = allocvector(nc1);
    itemp = 1;
    beta = allocvector(itemp);
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
    
    // test house function
    itemp = 0;
    house(X, v, beta, &nr, &itemp, &itemp);
    
    printf("Matrix X is \n");
    printm(X,&nr,&nc1);
    
    itemp=1;
    printf("The vector v for the first column j=1 is \n");
    printm(v,&itemp,&nc1);
    
    printf("beta = %.5f \n",beta[0]);
    puts("");
    
    printf("Calling the householder QR algorithm... \n");
    
    qr_household(X, &nr, &nc1);
    
    printf("The QR transformed matrix X is \n");
    printm(X,&nr,&nc1);
    
    printf("Call the backward accumulation algorithm to recover Q \n");
    back_accum(X, Q, &nr, &nc1);
    printm(Q,&nr,&nr);
    
    // fill in X and Y for least squares solution
    for (i=0; i<nr; i++) {
        for (j=0; j<(nc-1); j++) {
            X[j*nr+i]=dat[j*nr+i];
            Y[i]=dat[(nc-1)*nr+i];
        }
    }
    
    printf("Solving the least squares problem \n");
    linreg_hh(X, Y, &nr, &nc1);
    printv(Y, &nc1);
    
    
}

/***********************************************************************************/




// Reference: Algorithm 5.1.1 Matrix Computations, Golub & Van Loan 2013
// computes a household vector from a  long vector A and returns v and beta
// by subsetting the ia'th to m'th element of column ic in A
// A - mxn long vector, v - vector of length m, beta - scalar of length 1,

void house(double *A, double *v, double *beta, int *m, int *ia, int *ic)
{
    int i;
    double sig, mu, x1;
    v[*ia] = 1.0;
    x1 = A[(*ic)*(*m) + *ia]; //extract first element of x
    
    sig = 0.0; // initialise sig
    for (i=(*ia+1); i<*m; i++) { //construct initial nu
        sig+= A[(*ic)*(*m) + i]*A[*ic*(*m) + i];
        v[i] = A[(*ic)*(*m) + i];
    } // end for
    
    if (sig==0.0 & x1>=0) {
        beta[0] = 0.0;
    }else if(sig==0.0 & x1<0)
    {
        beta[0] = -2.0;
    }else{
        mu = sqrt((x1*x1) + sig);
        if (x1<=0.0) {
            v[*ia] = x1 -  mu;
        }else{
            v[*ia] = -sig/(x1 + mu);
        } // end if else
        beta[0] = 2*(v[*ia]*v[*ia])/(sig + (v[*ia]*v[*ia]));
        // compute v=v/v[1]
        for (i=*ia; i<*m; i++) {
            v[i] = v[i]/v[*ia];
        } // end for
    }// end if else if else
} // end void function

/***********************************************************************************/

// Reference: Algorithm 5.2.1 Matrix Computations, Golub & Van Loan 2013
void qr_household(double *A, int *m, int *n)
{
    int j,ir,ic,ik,lb=1;
    double *beta, *v, *vtemp, sum;
    
    // allocate memory
    beta = allocvector(lb);
    v = allocvector(*m);
    vtemp = allocvector(*m); // for temporary matrix multiplication values
    
    // main for loop
    for (j=0; j<*n; j++) {
        house(A, v, beta, m, &j, &j);
        
        ic=j; // index column
            for (ir=j; ir<*m; ir++) {
                sum = 0.0;
                for (ik=j; ik<*m; ik++) {
                    if (ik!=ir) {
                        sum+= (-beta[0]*v[ir]*v[ik])*A[ic*(*m)+ik];
                    }else{
                        sum+= (1-beta[0]*v[ir]*v[ik])*A[ic*(*m)+ik];
                    } // end if (ik!=ir) else{}
                } // end for ik
                vtemp[ir]=sum;
            } // end for ir
            for (ir=j; ir<*m; ir++) {
                A[ic*(*m) + ir] = vtemp[ir]; // fill in ic'th column of A
            }
        
        if (j<(*m-1)) {
            for (ir=(j+1); ir<*m; ir++) {
                A[j*(*m)+ir] = v[ir];//v[ir-j];
            }
        }// end if
        
    } // end for j
} // end void function

/***********************************************************************************/
// Reference: Algorithm (equation 5.1.5) Matrix Computations, Golub & Van Loan 2013
// algorithm for backward accumulation. this algorithm recovers the Q matrix from qr_household()

void back_accum(double *A, double *Q, int *m, int *n)
{
    int j ,i,ik,ir;
    double sig, beta, sum, *v, *vtemp;
    
    // allocate memory
    v = allocvector(*m);
    vtemp = allocvector(*m); // for temporary matrix multiplication values
    
    // initialise Q to identity matrix
    
    for (i=0; i<*m; i++) {
        for (j=0; j<=i; j++) {
            if (j!=i) {
                Q[j*(*m)+i] = 0.0;
                Q[i*(*m)+j] = 0.0; //by symmetry
            }else{
                Q[i*(*m)+j]=1.0;
            }// end if()
        }// end for j
    } //end for i
    
    
    for (j=(*n-1); j>=0; j--) {
        v[j] = 1.0;
        sig = 0.0;
        for (i=(j+1); i<*m; i++) {
            v[i] = A[j*(*m) + i];
            sig+= v[i]*v[i];
        }
        beta = 2/(1.0 + sig); // compute beta_j
        
        for (ir=j; ir<*m; ir++) {
            sum = 0.0;
            for (ik=j; ik<*m; ik++) {
                    sum+= (beta*v[ir]*v[ik])*Q[j*(*m)+ik];
            } // end for ik
            vtemp[ir]=sum;
        } // end for ir
        for (ir=j; ir<*m; ir++) {
            Q[j*(*m) + ir] = Q[j*(*m) + ir] - vtemp[ir]; // fill in ic'th column of A
        }
        
    } // end for j
}

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
// Reference: Algorithm (equation 5.3.2) Matrix Computations, Golub & Van Loan 2013
// algorithm for Householder LS Solution. This algorithm solves the Least Squares problem.

void linreg_hh(double *X, double *Y, int *nrx, int *ncx)
{
    int ncY = 1,i,j;
    double b, sig,zig, *v, tmp;
    
    // allocate memory
    v = allocvector(*nrx);
    
    qr_household(X, nrx, ncx); // QR factorisation; X overwritten with R in upper triangular
    for (j=1; j<*ncx; j++) {
        v[j]=1.0;
        for (i=(j+1); i<*nrx; i++) { // fill in v
            v[i] = X[j*(*nrx)+i];
        } // end for i
        sig=0.0; zig=0.0;
        for (i=j; i<*nrx; i++) { // compute beta
            sig += (v[i]*v[i]);
            zig += v[i]*Y[i];
        } // end for i
        b = (2/sig)*zig; // compute coefficient of v
        for (i=j; i<*nrx; i++) { // compute beta
            tmp = Y[i] - (b*v[i]);
            Y[i] = tmp;
        } // end for i
    } // end for j
    
    // solve for the coefficients
    backsolve(X, Y, ncx); // solution overwrites first ncx elements in Y
}


//
//  chol_decomp.c
//  
//
//  Created by Emmanuel Tsyawo on 4/25/19.
//

/* Compile using
 gcc -c chol_decomp.c matply.c utils.c
 gcc -o execCD chol_decomp.o matply.o utils.o
 ./execCD
 */

#include "matply.h"
#include "utils.h"
void chol_decomp(double *A, int *n);
int main()
{
    int n = 6,i,j;
    double *A2;
    
    // allocate memory
    A2 = allocvector(n*n);
    
    printf("\n");
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
        3.446,88.798,5.796,-2.039,-18.366,-15.693,
        2.686,5.796,85.911,7.203,-16.848,-16.081,
        6.798,-2.039,7.203,96.208,-1.853,-1.8,
        2.59,-18.366,-16.848,-1.853,90.32,7.953,
        -19.055,-15.693,-16.081,-1.8,7.953,109.126};
    
    // cholesky decomposition:
    printf("Calling the cholesky decomposition algorithm on matrix A = \n");
    printm(A,&n,&n);
    chol_decomp(A, &n);
    printf("The elements of matrix decomposed A are: \n" );
    printm(A,&n,&n);
    
    // fill in zeros
    
    for (i=0; i<(n-1); i++) {
        for (j=(i+1); j<n; j++) {
            A[j*(n)+i] = 0.0;
        }
    }
    printf("The elements of matrix decomposed A are: \n" );
    printm(A,&n,&n);
    
    printf("Recover matrix A\n");
    trans(A, A2, &n,&n);
    matply_sym(A2, A, &n, &n);
    printm(A,&n,&n);
}


// cholesky decomposition, lower triangular L is the decomposition
void chol_decomp(double *A, int *n)
{
    int j,l,i;
    double sig, *v;
    
    // allocate memory to v
    v = allocvector(*n);
    
    for (j=0; j<*n; j++) {
        if (j>0) {
            for(i=j;i<*n;i++) {
            sig = 0.0; //initialise sig
            for (l=0; l<=(j-1); l++) {
                sig+= A[l*(*n)+i]*A[l*(*n)+j];
            }// end for l
                v[i] = A[j*(*n)+i]-sig;
        } // end for i
            // pass v to A
            for (i=j; i<*n; i++) {
                A[j*(*n)+i]=v[i];
            } // end for i
        } // end if
        for (i=j; i<*n; i++) {
            A[j*(*n)+i] = A[j*(*n)+i]/sqrt(A[j*(*n)+j]);
        } // end for i
    } // end for j
}

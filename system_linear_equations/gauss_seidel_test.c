/* 
  Emmanuel S. Tsyawo 
  estsyawo@temple.edu,  estsyawo@gmail.com
  December 12, 2018
  Gauss Seidel algorithm for system of linear equations AX = b
  Reference: https://en.wikipedia.org/wiki/Gauss%E2%80%93Seidel_method 
*/
/* Compile using
 gcc -c gauss_seidel_test.c matply.c
 gcc -o execGauss_Seidel_test gauss_seidel_test.o matply.o
 ./execGauss_Seidel_test
 */

#include "matply.h"

// declare function prototypes.
void gauss_seidel( double *A, double *b, double *phi, double *dev, int *n);

// main function
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6;
    double *dev, *phi;
    dev = malloc(nra*sizeof(double));
    phi = malloc(nra*sizeof(double));

    printf("\n");
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
                          3.446,88.798,5.796,-2.039,-18.366,-15.693,
                          2.686,5.796,85.911,7.203,-16.848,-16.081,
                          6.798,-2.039,7.203,96.208,-1.853,-1.8,
                          2.59,-18.366,-16.848,-1.853,90.32,7.953,
                          -19.055,-15.693,-16.081,-1.8,7.953,109.126};

    double X[6] = {-1.0,-0.5,0.0,0.5,1.0,1.5};
    double b[6] = {-111.501,-90.771,-42.952,37.773,107.917,197.644};

    printf("The elements of matrix A are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", A[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");
    printf("The elements of vector b are: \n");
    for (i=0; i<nra; i++) {
        printf("b[%d] = %.3f \n",i,b[i]);
    }
    printf("\n");
    printf("The vector X of values to be solved for are: \n");
    for (i=0; i<nra; i++) {
        printf("X[%d] = %.1f \n",i,X[i]);
    }

    printf("Calling the Gauss-Seidel algorithm to solve for Z in AZ=b \n");

    gauss_seidel( A, b, phi, dev, &n);
    printf("\n");
    printf("The solution is the vector: \n");

    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,phi[i]);
    }


}


void gauss_seidel( double *A, double *b, double *phi, double *dev, int *n){
    double sig, vl, tol, mxdev;
    int i, j;
    tol=1e-10;
    for(;;) { //begin do while loop
        for(i=0; i<*n; i++){
            sig = 0.0;

            for(j=0; j<*n; j++ ){
                if(j!=i){
                    sig = (double) sig + A[j*(*n) + i]*phi[j];
                } // end if
            }// end for j

            if(absval(A[i*(1+(*n))])<tol){
                break;
                // add warning for non-dominant diagonal term
            }
            vl = phi[i];
            phi[i] = (double) ((b[i]-sig) - phi[i])/A[i*(1+(*n))];
            dev[i] = absval((phi[i]-vl));
        }// end for i
        mxdev = (double) max(dev,n);
        if(mxdev<=tol)
            break;
    }
}

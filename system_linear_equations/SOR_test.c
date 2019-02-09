/*
 Emmanuel S. Tsyawo
 estsyawo@temple.edu,  estsyawo@gmail.com
 December 12, 2018
 SOR for system of linear equations AX = b
 Reference: https://en.wikipedia.org/wiki/Successive_over-relaxation
 */
/* Compile using
 gcc -c SOR_test.c matply.c
 gcc -o execSOR_test SOR_test.o matply.o
 ./execSOR_test
 */

#include "matply.h"

// prototype function declaration
void SOR( double *A, double *b, double *phi, int *n);
// main function
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6;
    double *dev, *phi;
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
        phi[i]=0.0; //initialise phi to zeros
    }
    
    printf("Calling the SOR algorithm to solve for Z in AZ=b \n");
    
    SOR( A, b, phi, &n);
    printf("\n");
    printf("The solution is the vector: \n");
    
    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,phi[i]);
    }
    
    
}


// SOR void function
void SOR( double *A, double *b, double *phi, int *n){
    double sig, vl, tol, dev, mxdev, w, zi;
    int i, j, iter, maxiter;
    tol=1e-10;
    w = 1.25; // this value can be adjusted on the interval (0,2)
    maxiter = 1000;
    iter = 0;
    for(;;) { //begin do while loop
        iter +=1;
        for(i=0; i<*n; i++){
            sig = 0.0;
            
            for(j=0; j<*n; j++ ){
                if(j!=i){
                    sig = (double) sig + A[j*(*n) + i]*phi[j];
                } // end if
            }// end for j
            
            if(absval(A[i*(1+(*n))])<tol){
                break;
            }
            
            vl = (double) w*(((b[i]-sig)/A[i*(1+(*n))]) - phi[i]);
            zi = phi[i];
            phi[i] = phi[i] + vl;
            dev = absval((phi[i]-zi));
	    if(dev>mxdev){
	      mxdev = dev;
	     }
        }// end for i
        if(mxdev<=tol){
            break;
    }
        if (iter>=maxiter) {
            printf("Maximum number of iterations reached. Algorithm failed to converge.");
            break;
        }
        
    }
}

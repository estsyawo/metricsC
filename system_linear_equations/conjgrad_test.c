/*
! Emmanuel S. Tsyawo 
! estsyawo@temple.edu,  estsyawo@gmail.com
! January 1, 2019
! Conjugate Gradient algorithm for system of linear equations AX = b
! Reference: https://en.wikipedia.org/wiki/Conjugate_gradient_method
*/
/* Compile using
 gcc -c conjgrad_test.c matply.c
 gcc -o execConjGrad_test conjgrad_test.o matply.o
 ./execConjGrad_test
 */
#include "matply.h"
// conjugate gradient method for solving a system of linear equations AX=b

// prototype function declaration

void conjgrad(double *A, double *b, double *X, int *n);

// main program
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6;
    double *Z;

    Z = malloc(nra*sizeof(double)); // vector to store the solution

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

    printf("\n");
    printf("Calling the Conjugate gradient algorithm to solve for Z in AZ=b \n");

    conjgrad(A, b, Z, &n);
    printf("\n");
    printf("The solution is the vector: \n");

    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,Z[i]);
    }

}

// function for the conjugate gradient algorithm
void conjgrad(double *A, double *b, double *X, int *n)
{
    double *r, *p, *vec;
    double alf, beta, tol, dev, dpr1, dpr0;
    int ncX, i,k,maxiter;

    // set parameters
    ncX = 1; //number of columns in X
    tol = 1e-7;
    maxiter = 1000; // maximum number of iterations allowed.

    // allocate vectors r, p and vector holder vec
    r = malloc((*n)*sizeof(double));
    p = malloc((*n)*sizeof(double));
    vec = malloc((*n)*sizeof(double));

    k = 0;
    for(i=0;i<*n;i++){ // fill in initial r
        X[i] = 0.0; // initialise X to zero => A*X = 0 => initial r=b
        r[i] = b[i];
        p[i] = r[i]; //initialise p with r
        //printf("vec[%d]=%.3f, r[%d]=%.3f\n",i,vec[i],i,r[i]);
    }//end for i

    dpr0=dotprod(r, r, n); // take dot product
    // initialise main while loop
    for(;;){

        matply(A, p, vec, n, n, &ncX); // take matrix product vec=A*p
        dpr1=dotprod(p, vec, n);// take dot product, store in dpr1 temporarilly

        alf = dpr0/dpr1;
        //printf("alf = %.3f, k = %d\n",alf,k);

        //update X
        for(i=0;i<*n;i++){ // update X, r
            X[i] += alf*p[i];
            r[i] += -alf*vec[i] ; // recall vec=Ap currently
            vec[i] = absval(r[i]); // pass |r| to vec
        }

        dev = max(vec, n); //compute infinity norm of r
        //printf("dev = %.3f, k = %d\n",dev,k);

        // checking for convergence
        if(dev<tol){
            break;
        }
        if(k>=maxiter){
            printf("Maximum number of iterations reached. Conjugate gradient algorithm failed to converge.\n");
            break;
        }

        dpr1 = dotprod(r,r,n); // take dot product dpr1
        beta = dpr1/dpr0;

        for(i=0;i<*n;i++){// update p
            p[i] = r[i] + beta*p[i];
        }

        k += 1; //update k
        dpr0=dpr1; // update drp0

    }// end for(;;)

    printf("Number of iterations = %d\n",k);
    // free allocated memory
    free(r);
    free(p);
    free(vec);
}

#include "matply.h"

// conjugate gradient method for solving a system of linear equations AX=b

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
    }//end for i

    dpr0=dotprod(r, r, n); // take dot product
    // initialise main while loop
    for(;;){

        matply(A, p, vec, n, n, &ncX); // take matrix product vec=A*p
        dpr1=dotprod(p, vec, n);// take dot product, store in dpr1 temporarily

        alf = dpr0/dpr1;

        //update X
        for(i=0;i<*n;i++){ // update X, r
            X[i] += alf*p[i];
            r[i] += -alf*vec[i] ; // recall vec=Ap currently
            vec[i] = absval(r[i]); // pass |r| to vec
        }

        dev = max(vec, n); //compute infinity norm of r

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

        k += 1; //update number of iterations k
        dpr0=dpr1; // update drp0

    }// end for(;;)

    // free allocated memory
    free(r);
    free(p);
    free(vec);
}

/*Compile using the following files: conjgrad.c matply.c*/

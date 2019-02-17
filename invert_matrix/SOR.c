#include "matply.h"

void SOR( double *A, double *b, double *phi, int *n){
    double sig, vl, tol, dev, mxdev, w;
    int i, j, iter, maxiter;
    tol=1e-7;
    w = 1.0; // this value can be adjusted on the interval (0,2)
    maxiter = 100;
    iter = 0;
    for(;;) { //begin do while loop
        iter +=1;
        for(i=0; i<*n; i++){
            sig = 0.0;
            dev = 0.0; mxdev = 0.0; // reset in order to check convergence
            for(j=0; j<*n; j++ ){
                if(j!=i){
                    sig = (double) sig + A[j*(*n) + i]*phi[j];
                } // end if
            }// end for j
            
            if(absval(A[i*(1+(*n))])<tol){
                break;
            }
            
            vl = (double) w*(((b[i]-sig)/A[i*(1+(*n))]) - phi[i]);
            phi[i] = phi[i] + vl;
            dev = absval(vl);
            if(dev>mxdev){
                mxdev = dev; // update infinity norm of deviations
            }
        }// end for i
        if(mxdev<=tol){ // check for convergence
            break;
        }
        if (iter>=maxiter) {
            printf("Warning: Maximum number of iterations reached.\n");
            break;
        }
        
    }
}
/*
 Compile using: SOR.c matply.c
 */

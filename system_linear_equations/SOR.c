#include "matply.h"

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
            printf("Maximum number of iterations reached. Algorithm failed to converge.\n");
            break;
        }
        
    }
}
/*
 Compile using: SOR.c matply.c
 */

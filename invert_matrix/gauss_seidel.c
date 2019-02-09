#include "matply.h"

void gauss_seidel( double *A, double *b, double *x, int *n){
	double sig, vl, tol, dev, mxdev;
	int i, j, iter, maxiter;
	tol=1e-10;
	iter = 0;
	maxiter = 1000;

for(;;) { //begin do while loop
	iter +=1;
	for(i=0; i<*n; i++){
	  sig = 0.0; dev=0.0; mxdev=0.0;
		for(j=0; j<*n; j++ ){
		if(j!=i){
		sig = (double) sig + A[j*(*n) + i]*x[j];
		    } // end if
	}// end for j

	if(absval(A[i*(1+(*n))])<tol){
		printf("Algorithm stopped: non-dominant diagonal term");
	break;
	}
	vl = x[i];
	x[i] = (double) ((b[i]-sig) - x[i])/A[i*(1+(*n))];
	dev = absval((x[i]-vl));
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

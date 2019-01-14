#include "matply.h"

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
// checking convergence
//return *phi;
}

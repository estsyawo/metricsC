//
//  matrix_decompositions.c
//  
//
//  Created by Emmanuel Tsyawo on 3/11/19.
// Reference: https://en.wikipedia.org/wiki/QR_decomposition

#include "matply.h"
#include "utils.h"
// QR Decomposition by householder reflections
// Dimensions: A - mxn, Q - mxm, Qk - mxm
void qr_hh(double *A, double *Q, double *Qk, int *m, int *n)
{
    int i,j,t,k,l, tm;
    double alf, un, *u, sig;
    // allocate memory
    u = allocvector(*m);
    // set parameters
    tm = min_2(*m-1,*n); // total number of loops
    
    
    // begin loops
    for (t=0; t<tm; t++) {
        sig = 0.0;
        for (l=t; l<*m; l++) {
            sig += (A[t*(*m)+l]*A[t*(*m)+l]);
        }
        alf = sqrt(sig); // compute the norm of x
        // compute u<-- x-alf*e
        u[t] = A[t*((*m)+1)]-alf; //compute first element of u
        
        for (k=t+1; k<*m; k++) {
            u[k]=A[t*(*m)+k]; //fill in the rest of u
        }
        
        sig = 0.0;
        for (l=t; l<*m; l++) {
            sig += (u[l]*u[l]);
        }
        un = sqrt(sig); // un<-- ||u||
        
        if (un < 1e-10) { // check for rank deficiency
            printf("Error: Floating point error; division by zero. Possible rank deficiency. \n");
            break;
        }
        
        //compute v<--u/||u||
        for (k=t; k<*m; k++) {
            u[k]=u[k]/un; //compute v<--u/||u||
        }
        
        // Compute Qk and fill in 1's and zeros in left upper portion
        for (i=t; i<*m; i++) {
            for (j=t; j<i; j++) {
                Qk[j*(*m)+i]=-2*u[i]*u[j];
                Qk[i*(*m)+j]=Qk[j*(*m)+i]; // by symmetry
                if (t>0) {
                    Qk[(t-1)*(*m)+i-1]=0.0; // fill in zeros for off-up-diag elements
                    Qk[(i-1)*(*m)+(t-1)]=0.0;
                }
            } // end for j
            Qk[i*(*m+1)]=1-(2*(u[i]*u[i]));
        } // end for i
        

        
        if (t>0) {
            Qk[(t-1)*(*m)+ *m-1]=0.0; // fill in zeros for off-up-diag elements
            Qk[(*m-1)*(*m)+(t-1)]=0.0;
            Qk[(t-1)*(*m+1)]=1.0; // fill in 1 for upper diag element
        }
        printf("iter = %d ",t); printf("Matrix Qk is "); printm(Qk,m,m);
        matply(Qk, A, A, m, m, n); //update R i.e A<-- Qk*A; A destroyed to store R
        printf("iter = %d ",t); printf("Matrix R is "); printm(A,m,n);
        
        if (t>0) {
            matply_xyt(Q, Qk, Q, m, m, m); //update Q<-- Qt-1'Qt'
        }else{
            trans(Qk, Q,m,m); // Q<-- Q'
        }
    }// end for (t=0; t<tm; t++)
}

/*
 Emmanuel S. Tsyawo
 estsyawo@temple.edu,  estsyawo@gmail.com
 February 08, 2019
 A wrapper function for solvers of system of linear equations
 */
/* Compile using
 gcc -c SOR_test.c matply.c
 gcc -o execSOR_test SOR_test.o matply.o
 ./execSOR_test
 */
//**************************************************************************************//
// This file pools together all solvers for systems of linear equations
//**************************************************************************************//

#include "solve.h"
#include "matply.h"

/*  
  Input:
	A - nxn matrix
	b - nx1 vector
	X - nx1 vector of unknowns to be solved for
	type - a character for type of solver. solvers currently available are
	"conjgrad", "gauss_seidel", and "SOR"
*/
//**************************************************************************************//
void solve(double *A, double *b, double *X, int *n, char *type)
{
  switch(*type){
    case 'C':
	  conjgrad(A, b, X, n);
	  break;
	case 'G':
	  gauss_seidel(A, b, X, n);
	  break;
	case 'S':
	  SOR(A, b, X, n);
	  break;
    case 'Q':
      lsqqr( A, b, X, n, n);
      break;
	default :
	  printf("Invalid type of solver.\n");
  }	

}
//**************************************************************************************//

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


//**************************************************************************************//

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

//**************************************************************************************//

void SOR( double *A, double *b, double *phi, int *n){
    double sig, vl, tol, dev, mxdev, w, zi;
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

//**************************************************************************************//


/* QR algorithm accessed from:
 Source: http://www.math.umd.edu/~mariakc/teaching-2/householder.c
 on February 22, 2019
 */

void lsqqr( double *a, double *b, double *coef, int *nrow,  int *ncol) {
    int i,m,l;
    double normx, aaux, baux;
    double aux[*ncol], u[*nrow];
    double c[*nrow]; /* weight matrix */
    
    /* put larger weights at the ends interval to reduce the error at the ends of the interval */
    /* Ax = b   <==> CAx = Cb where C = diag{c[0],c[1],...,c[nrow]}  */
    for( m=0; m<*nrow; m++ ) c[m]=1.0;
    for( m=0; m<min_2((*nrow)/2,10); m++ ) c[m]=(10-m);
    for( m=(*nrow)-1; m>max_2((*nrow)/2,(*nrow)-9); m-- ) c[m]=(9-(*nrow)+m);
    /* compute CA and Cb */
    for( m=0; m<*nrow; m++ ) {
        for( i=0; i<*ncol; i++ ) {
            (*(a+i+m*(*ncol)))*=c[m];   /* a_{m,i} = a_{m,i} * c_m */
        }
        b[m]*=c[m];   /* b_m = b_m * c_m */
    }
    
    /* Start the QR algorithm via Householder reflections */
    for( i=0; i<*ncol; i++ ) {
        /*form vector u=House(a(i:nrow,i))*/
        normx=0.0;
        for( m=0; m<*nrow-i; m++ ) {
            aaux=*(a+(i+m)*(*ncol)+i);  /* aaux = a_{m+i,i} */
            normx+=aaux*aaux;  /* normx = normx + aaux*aaux */
        }
        aaux=*(a+i*(*ncol)+i);    /* aaux = a_{i,i} */
        u[0]=aaux+sgn(aaux)*sqrt(normx);
        for( m=1; m<*nrow-i; m++ ) u[m]=*(a+(i+m)*(*ncol)+i);  /* u[m] = a_{i+m,i} */
        normx=0.0;
        for( m=0; m<*nrow-i; m++ ) normx+=u[m]*u[m];   /* normx = normx + u[m]*u[m] */
        /* compute (I-2uu^t)a(i:nrow,i:ncol) */
        for( l=0; l<*ncol-i; l++ ) {
            aux[l]=0.0;
            for( m=0; m<*nrow-i; m++ )  {
                aux[l]+=u[m]*(*(a+(i+m)*(*ncol)+i+l));  /* aux[l] = aux[l] + u[m]*a_{i+m,i+l} */
            }
            aux[l]*=(-2.0/normx);     /* aux[l] = -2*aux[l]/normx */
            for( m=0; m<*nrow-i; m++ )  {
                *(a+(i+m)*(*ncol)+i+l)+=u[m]*aux[l];      /* a_{i+m,i+l} = a_{i+m,i+l} + u[m]*aux[l] */
            }
        }
        /* compute (Q^T)*b */
        baux=0.0;
        for( m=0; m<*nrow-i; m++ ) baux+=u[m]*b[i+m];  /* baux = baux + u[m]*b[i+m] */
        
        baux*=(-2.0/normx);
        for( m=0; m<*nrow-i; m++ ) b[i+m]+=baux*u[m];  /* b[i+m] = b[i+m] + baux*u[m] */
    }
    /* find the coefficients*/
    for( i=(*ncol)-1; i>=0; i-- ) {
        baux=0.0;
        for( m=i+1; m<*ncol; m++ ) {
            baux+=(*(a+i*(*ncol)+m))*coef[m];    /* baux = baux + a_{j,m}*coef[m] */
        }
        coef[i]=(b[i]-baux)/(*(a+i*(*ncol)+i));
    }
}




//**************************************************************************************//

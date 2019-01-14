// multiply two matrices
/*
This file contains functions for matrix and vector operations and other useful
functions
*/

# include <math.h>

/* xa - vectorised matrix A, 
  xb - vectorised matrix B,
  xab - product AB (vectorised)
  nra - number of rows in A
  nca - number of columns in B
  ncb - number of columns in B
*/


void matply(double *xa, double *xb,double *xab, int *nra, int *nca,int *ncb){
    double sum ;
    for(int i=0; i< *nra;i++){
        for (int j=0; j<*ncb; j++) {
            sum = 0.0;
            for (int k=0; k<*nca; k++) {
                sum =  (double) (sum + (xa[ k*(*nra)+i])*(xb[ j*(*nca)+k])) ;
            }
         xab[ j*(*nra)+i] = sum ;
        }
    }
}

/*
Take the transpose of a matrix a, store in atrans,
*/
void trans(double *a, double *atrans,int *nra,int *nca)
{
  int i,j;
    for(i=0;i<*nra;i++){
      for(j=0;j<*nca;j++){
	  atrans[i*(*nca) +j] = a[j*(*nra)+i];
	}
    }
}

// take the dot product two vectors a, b with length n each.
double dotprod(double *a, double *b, int *n)
{
    int i;
    double ans;
    ans = 0.0;
    for (i=0; i<*n; i++) {
        ans += a[i]*b[i];
    }
    return ans;
}

// find the maximum in a vector x of length n
double max(double *x, int *n)
{
    double xmax;
    for(int i=1; i<*n; i++){
        if(x[i]>x[i-1])
            xmax=x[i];
        else
            xmax=x[i-1];
    }
    return xmax;
}

// find the minimum in a vector x of length n
double min(double *x, int *n)
{
    double xmin;
    for(int i=1; i<*n; i++){
        if(x[i]<x[i-1])
            xmin=x[i];
        else
            xmin=x[i-1];
    }
    return xmin;
}


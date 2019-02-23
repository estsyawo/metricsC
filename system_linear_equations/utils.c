//  utils.c
//  
//
//  Created by Emmanuel Tsyawo on 2/9/19.
//
#include "utils.h"
//*************************************************************************************//
// print output as vector or matrix

// print as vector
void printv(double *vec, int *n)
{
    int i;
    double v;
    puts(" ");
    for (i=0; i<*n; i++) {
        v = vec[i];
        if (v>=0.0) {
            printf(" %.2f \n",v);
        }else{
            printf("%.2f \n",v);
        }
    }
    puts(" ");
}

// print as matrix for doubles
void printm(double *mat, int *nr, int *nc)
{
    int i,j;
    double v;
    puts(" ");
    for (i=0; i<*nr; i++) {
        for (j=0; j<*nc; j++) {
            v = mat[j*(*nr) + i];
            if (v>=0.0) {
                printf(" %.2f ",v);
            }else{
                printf("%.2f ",v);
            }
        }
        puts(" ");
    }
    puts(" ");
}
// print as matrix for integers
void printm_int(int *mat, int *nr, int *nc)
{
    int i,j,v;
    puts(" ");
    for (i=0; i<*nr; i++) {
        for (j=0; j<*nc; j++) {
            v = mat[j*(*nr) + i];
            if (v>=0) {
                printf(" %d ",v);
            }else{
                printf("%d ",v);
            }
        }
        puts(" ");
    }
    puts(" ");
}

//*************************************************************************************//
// allocate a vector of size n to a pointer
double *allocvector(int n)
{
    double *vec;
    vec=(double *) malloc((size_t)((n)*sizeof(double)));
    if (!vec) {
        printf("Failure to allocate vector\n");
        exit(1);
    }
    return vec;
}
//*************************************************************************************//

// initialise vector vec of length n to all scalar a's
void init_vec (double *vec, int *n, double *a)
{
    int i;
    for (i=0; i<*n; i++) {
        vec[i] = *a;
    }
}

//*************************************************************************************//
// compare functions for use in qsort()

// compare function for an array of integers
int cmpfunc_int(const void * a, const void * b) {
    return ( *(int*)a - *(int*)b );
}

// compare function for an array of doubles
int cmpfunc_d (const void * a, const void * b)
{
    if (*(double*)a > *(double*)b)
        return 1;
    else if (*(double*)a < *(double*)b)
        return -1;
    else
        return 0;
}

// a wrapper for qsort() for double
void sort(double *x, int n)
{
    qsort(x, n, sizeof(double), cmpfunc_d);
}

// a wrapper for qsort() for integers
void sort_int(int *x, int n)
{
    qsort(x, n, sizeof(int), cmpfunc_int);
}
//*************************************************************************************//
/*
 Probability density and cumulative density functions
 d_ denotes the pdf and p_ denotes the cdf
 */

// logit
double d_logit(double x)
{
    return exp(x)/pow((exp(x)+1.0),2);
}

double p_logit(double x)
{
    return 1.0/(1.0+exp(-x));
}

// normal density
// mu - mean, sig - standard deviation
double d_norm(double x, double mu, double sig)
{
    double pi2 = 8.0 * atan(1.0);
    double ans = (1.0 / (sig * sqrt(pi2))) *
    exp(-(x - mu)*(x - mu) / (2.0 * sig * sig));
    return ans;
}

double p_norm(double x, double mu, double sig)
{
    return 0.5*erfc(-M_SQRT1_2*(x-mu)/sig);
}

//*************************************************************************************//

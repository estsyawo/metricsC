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
            printf(" %3.3f \n",v);
        }else{
            printf("%3.3f \n",v);
        }
    }
    puts(" ");
}

// print as matrix
void printm(double *mat, int *nr, int *nc)
{
    int i,j;
    double v;
    puts(" ");
    for (i=0; i<*nr; i++) {
        for (j=0; j<*nc; j++) {
            v = mat[j*(*nr) + i];
            if (v>=0.0) {
                printf(" %3.3f ",v);
            }else{
                printf("%3.3f ",v);
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
    vec=(double *) malloc((size_t)((n+1)*sizeof(double)));
    if (!vec) {
        printf("Failure to allocate vector\n");
        exit(1);
    }
    return vec;
}
//*************************************************************************************//

// initialise vector vec of length n to all a's
void init_vec (double *vec, int *n, double *a)
{
    int i;
    for (i=0; i<*n; i++) {
        vec[i] = *a;
    }
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

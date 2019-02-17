//
//  utils.h
//  
//
//  Created by Emmanuel Tsyawo on 2/9/19.
//
// This function contains motley fundamental functions that are called into several routines.

#include <stdio.h>
#include <math.h>
#include <stdlib.h>

// print matrix/vectors to terminal
void printv(double *vec, int *n);
void printm(double *mat, int *nr, int *nc);


// probability densities and distributions
double d_logit(double x);
double p_logit(double x);
double d_norm(double x, double mu, double sig);
double p_norm(double x, double mu, double sig);

// memory allocation
double *allocvector(int n);

// vector/matrix initialisation
void init_vec (double *vec, int *n, double *a);

// for reading .txt data
double *read_txt(char *datname, int *nr, int *nc );
void dat_read_prep(char *datname, double *X, double *Y, int *nr, int *nc );

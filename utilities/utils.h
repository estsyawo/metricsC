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

double *allocvector(int n);
double *read_txt(char *datname, int *nr, int *nc );
double d_logit(double x);
double p_logit(double x);
double d_norm(double x, double mu, double sig);
double p_norm(double x, double mu, double sig);

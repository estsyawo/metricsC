//
//  score_hessian.h
//  
//
//  Created by Emmanuel Tsyawo on 2/11/19.
//

void scorehess_logit(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat);
void new_raph_logit(double *beta, double *y, double *x, int *nrx, int *ncx);
void scorehess_probit(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat);
void new_raph_probit(double *beta, double *y, double *x, int *nrx, int *ncx);
void scorehess_poisson(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat);
void new_raph_poisson(double *beta, double *y, double *x, int *nrx, int *ncx);

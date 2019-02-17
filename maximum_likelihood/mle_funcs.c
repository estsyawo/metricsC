//
//  score_hessian.c
//  
//
//  Created by Emmanuel Tsyawo on 2/11/19.
//
/* Purpose: This functions contain likelihood functions (like_), score/gradient functions
 (score_), and hessian functions (hess_) for estimation and inference via maximum likelihood
 scorehess_ void functions compute both the score and the hessian
*/

#include "matply.h"
#include "utils.h"
#include "solve.h"

//************************************************************************************************//
// The logit model
// logit score and hessian functions

void scorehess_logit(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat)
{
    double pv, lc, dv, tmp;
    int i,j, itp,nlh;
    
    pv=0.0; nlh = (*ncx)*(*ncx);
    init_vec(score, ncx, &pv); // initialise the score vector to zeros
    init_vec(hessian, &nlh, &pv); // initialise the hessian matrix to zeros
    
    itp = 1; // set as number of rows in a row vector x[i,]
    for (i=0; i<*nrx; i++) {
        for (j=0; j<*ncx; j++) {
            tempvec[j] = x[j*(*nrx) + i]; // extract x[i,] into tempvec
        }
        matply_sym(tempvec, tempmat, &itp, ncx); // compute x'x and store in tempmat
        lc = dotprod(tempvec, beta, ncx); // take linear combination x[i,]beta
        pv = p_logit(lc); // logit cdf of lc
        dv = d_logit(lc); // logit pdf of lc
        dv = -1.0*dv; // negative pdf for the hessian
        pv = y[i] - pv; // for the score
        matply_ax(tempvec, &pv, tempvec, ncx); // compute score for observation i and store in tempvec
        vecadd(tempvec, score, score, ncx); // update score with observation i's score
        
        matply_ax(tempmat, &dv, tempmat, &nlh); // compute hessian for i
        vecadd(tempmat, hessian, hessian, &nlh); // update hessian
    }
}

// solving the logit model using the Newton-Raphson algorithm
// beta - starting values, y , x - design matrix (include 1's for intercept),
// nrx - number of observations, ncx - length of beta

void new_raph_logit(double *beta, double *y, double *x, int *nrx, int *ncx)
{
    double tol, dev, *tempvec, *tempmat, *score, *hessian ;
    int k,maxk;
    
    tol = 1e-7;
    k = 0, maxk = 20;
    
    // allocate memory for
    tempvec=allocvector(*ncx); // allocate memory for a vector of length ncx
    tempmat=allocvector((*ncx)*(*ncx)); // allocate memory
    score = allocvector(*ncx); // allocate memory to score
    hessian = allocvector((*ncx)*(*ncx)); // allocate memory to hessian
    
    for(;;){
        k+=1; // update iteration
    scorehess_logit(score, hessian, beta, y, x, nrx, ncx, tempvec, tempmat); // compute score and hessian
    // solve change in X
    conjgrad(hessian, score, tempvec, ncx); // use conjugate gradient solver, difference in tempvec
    dev = norm_max(tempvec,ncx); // compute norm of Delta beta
        if (dev<=tol) {
            break;
        }
    vecsub(beta, tempvec, beta, ncx);
        if (k>=maxk) {
            printf("Warning: maximum number of iterations reached.");
            break;
        }
        printf("Iteration = %d ended. \n",k);
    }
    
    free(tempvec);
    free(tempmat);
    free(score);
    free(hessian);
}
//************************************************************************************************//


//************************************************************************************************//
// The probit model
// probit score and hessian functions

void scorehess_probit(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat)
{
    double pv, lc, dv, cvar, tmp,mu,sig;
    int i,j, itp,nlh;
    
    pv=0.0; mu = 0.0; sig = 1.0;
    nlh = (*ncx)*(*ncx);
    init_vec(score, ncx, &pv); // initialise the score vector to zeros
    init_vec(hessian, &nlh, &pv); // initialise the hessian matrix to zeros
    
    itp = 1; // set as number of rows in a row vector x[i,]
    for (i=0; i<*nrx; i++) {
        for (j=0; j<*ncx; j++) {
            tempvec[j] = x[j*(*nrx) + i]; // extract x[i,] into tempvec
        }
        matply_sym(tempvec, tempmat, &itp, ncx); // compute x'x and store in tempmat
        lc = dotprod(tempvec, beta, ncx); // take linear combination x[i,]beta
        pv = p_norm(lc,mu,sig); // logit cdf of lc
        cvar = pv*(1-pv);
        dv = d_norm(lc,mu,sig); // logit pdf of lc
        pv = dv*(y[i] - pv)/cvar; // for the score
        matply_ax(tempvec, &pv, tempvec, ncx); // compute score for observation i and store in tempvec
        vecadd(tempvec, score, score, ncx); // update score with observation i's score
        
        dv = -1.0*dv*dv/cvar; // negative pdf for the hessian
        matply_ax(tempmat, &dv, tempmat, &nlh); // compute hessian for i
        vecadd(tempmat, hessian, hessian, &nlh); // update hessian
    }
}

// solving the logit model using the Newton-Raphson algorithm
// beta - starting values, y , x - design matrix (include 1's for intercept),
// nrx - number of observations, ncx - length of beta

void new_raph_probit(double *beta, double *y, double *x, int *nrx, int *ncx)
{
    double tol, dev, *tempvec, *tempmat, *score, *hessian ;
    int k,maxk;
    
    tol = 1e-7;
    k = 0, maxk = 20;
    
    // allocate memory for
    tempvec=allocvector(*ncx); // allocate memory for a vector of length ncx
    tempmat=allocvector((*ncx)*(*ncx)); // allocate memory
    score = allocvector(*ncx); // allocate memory to score
    hessian = allocvector((*ncx)*(*ncx)); // allocate memory to hessian
    
    for(;;){
        k+=1; // update iteration
        scorehess_probit(score, hessian, beta, y, x, nrx, ncx, tempvec, tempmat); // compute score and hessian
        // solve change in X
        conjgrad(hessian, score, tempvec, ncx); // use conjugate gradient solver, difference in tempvec
        dev = norm_max(tempvec,ncx); // compute norm of Delta beta
        if (dev<=tol) {
            break;
        }
        vecsub(beta, tempvec, beta, ncx);
        if (k>=maxk) {
            printf("Warning: maximum number of iterations reached.");
            break;
        }
        printf("Iteration = %d ended. \n",k);
    }
    
    free(tempvec);
    free(tempmat);
    free(score);
    free(hessian);
}
//************************************************************************************************//



//************************************************************************************************//
// The poisson model

void scorehess_poisson(double *score, double *hessian, double *beta, double *y, double *x, int *nrx, int *ncx,double *tempvec, double *tempmat)
{
    double pv, lc, dv, tmp;
    int i,j, itp,nlh;
    
    pv=0.0; nlh = (*ncx)*(*ncx);
    init_vec(score, ncx, &pv); // initialise the score vector to zeros
    init_vec(hessian, &nlh, &pv); // initialise the hessian matrix to zeros
    //printv(score,ncx); printm(hessian,ncx,ncx);
    itp = 1; // set as number of rows in a row vector x[i,]
    for (i=0; i<*nrx; i++) {
        for (j=0; j<*ncx; j++) {
            tempvec[j] = x[j*(*nrx) + i]; // extract x[i,] into tempvec
        }
        // printv(tempvec,ncx);
        matply_sym(tempvec, tempmat, &itp, ncx); // compute x'x and store in tempmat
        lc = dotprod(tempvec, beta, ncx); // take linear combination x[i,]beta
        pv = exp(lc); // exp() of lc
        dv = -1.0*pv; // negative pdf for the hessian
        //dv = pv;
        pv = y[i] - pv; // for the score
        //printf("obs = %d, %3.3f\n",i,pv);
        matply_ax(tempvec, &pv, tempvec, ncx); // compute score for observation i and store in tempvec
        //printv(tempvec,ncx);
        vecadd(tempvec, score, score, ncx); // update score with observation i's score
        
        matply_ax(tempmat, &dv, tempmat, &nlh); // compute hessian for i
        vecadd(tempmat, hessian, hessian, &nlh); // update hessian
    }
}


void new_raph_poisson(double *beta, double *y, double *x, int *nrx, int *ncx)
{
    double tol, dev, *tempvec, *tempmat, *score, *hessian ;
    int k,maxk;
    
    tol = 1e-7;
    k = 0, maxk = 20;
    
    // allocate memory for
    tempvec=allocvector(*ncx); // allocate memory for a vector of length ncx
    tempmat=allocvector((*ncx)*(*ncx)); // allocate memory
    score = allocvector(*ncx); // allocate memory to score
    hessian = allocvector((*ncx)*(*ncx)); // allocate memory to hessian
    
    for(;;){
        k+=1; // update iteration
        scorehess_poisson(score, hessian, beta, y, x, nrx, ncx, tempvec, tempmat); // compute score and hessian
        conjgrad(hessian, score, tempvec, ncx);
        //printv(tempvec,ncx); printm(tempmat,ncx,ncx);
        vecsub(beta, tempvec, beta, ncx);
        dev = norm_max(tempvec,ncx); // compute norm of Delta beta
        if (dev<=tol) {
            break;
        }
        
        if (k>=maxk) {
            printf("Warning: maximum number of iterations reached.\n");
            break;
        }
        printf("Iteration = %d ended. \n",k);
    }
    
    free(tempvec);
    free(tempmat);
    free(score);
    free(hessian);
}

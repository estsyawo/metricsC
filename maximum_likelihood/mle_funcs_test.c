//
//  mle_funcs_test.c
//  
//
//  Created by Emmanuel Tsyawo on 2/12/19.
//
#include "mle_funcs.h"
#include "utils.h"


/* Compile using:
    gcc -c mle_funcs_test.c mle_funcs.c utils.c matply.c read_txt.c solve.c
    gcc -o execMle_F mle_funcs_test.o mle_funcs.o utils.o matply.o read_txt.o solve.o
    ./execMle_F
 */

int main()
{
    int nr=1000, nc = 6; // dimension of data set
    double *X, *Y, *score, *beta, *hessian, *tempvec,*tempmat,v;
    char *datname;
    
    double tbet[6] = {1.0648,0.5346,0.941,0.058,1.612,-0.5433}; // "True parameters" estimated in R
    Y = allocvector(nr); // allocate memory to Y
    tempvec=allocvector(nc); // allocate memory for a vector of length ncx
    tempmat=allocvector(nc*nc); // allocate memory
    score = allocvector(nc); // allocate memory to score
    beta = allocvector(nc); // allocate memory to starting values beta
    hessian = allocvector(nc*nc); // allocate memory to hessian
    
    //***********************************************************************//
    
    printf("Estimating the logit model \n");
    puts("");
    
    datname= "dat_binrep.txt"; // name of data set in folder
    // read in data
    X=read_txt(datname, &nr, &nc );
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nr, &nc );

    printf("Printing the parameter estimates to be solved for \n");
    printv(tbet,&nc);
    
    printf("Using the Newton-Raphson to solve the logit \n");
    
    // initialise beta to zeros
    v = 0.0;
    init_vec(beta, &nc, &v);
    new_raph_logit(beta, Y, X, &nr, &nc);
    printf("The solution is beta = \n");
    printv(beta,&nc);
    
    scorehess_logit(score, hessian, beta, Y, X, &nr, &nc,tempvec, tempmat);
    printf("Printing the score function \n");
    printv(score,&nc);
    printf("Printing the hessian function \n");
    printm(hessian,&nc,&nc);
    
    puts("");

    
    //***********************************************************************//
 
    printf("Estimating the probit model \n");
    puts("");
    // parameter vector of parameters estimated in R
     double tbet_probit[6] = {0.58085540, 0.26875450, 0.47809914, 0.03478482, 0.82715094, -0.26431253};
     // parameter vector of parameters estimated in R
     datname= "dat_binrep.txt"; // name of data set in folder
     // read in data
     X=read_txt(datname, &nr, &nc );
     
     // prepare data for regression
     dat_read_prep(datname, X, Y, &nr, &nc );
     
     printf("Printing the parameter estimates to be solved for \n");
     printv(tbet_probit,&nc);
     
     printf("Using the Newton-Raphson to solve the probit \n");
     
     // initialise beta to zeros
     v = 0.0;
     init_vec(beta, &nc, &v);
     new_raph_probit(beta, Y, X, &nr, &nc);
     printf("The solution is beta = \n");
     printv(beta,&nc);
     
     scorehess_probit(score, hessian, beta, Y, X, &nr, &nc,tempvec, tempmat);
     printf("Printing the score function \n");
     printv(score,&nc);
     printf("Printing the hessian function \n");
     printm(hessian,&nc,&nc);
     
     puts("");

    //***********************************************************************//

    printf("Estimating the poisson model \n");
    puts("");
    double tbet_pois[6] = {0.12106086,0.07241156,0.11947498,-0.01388003,0.17791073,-0.05883644}; // True parameters generating the data set
    
    datname= "dat_pois.txt"; // name of data set in folder
    // read in data
    X=read_txt(datname, &nr, &nc );
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nr, &nc );
    
    printf("Printing the parameter estimates to be solved for \n");
    printv(tbet_pois,&nc);
    
    printf("Using the Newton-Raphson to solve the poisson model \n");
    // initialise beta to zeros
    v = 0.0;
    init_vec(beta, &nc, &v);
    
    new_raph_poisson(beta, Y, X, &nr, &nc);
    printf("The solution is beta = \n");
    printv(beta,&nc);
    
    scorehess_poisson(score, hessian, beta, Y, X, &nr, &nc,tempvec, tempmat);
    printf("Printing the score function \n");
    printv(score,&nc);
    printf("Printing the hessian function \n");
    printm(hessian,&nc,&nc);
    puts("");
}


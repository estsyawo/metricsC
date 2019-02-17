//  linreg_test.c
//
//  Created by Emmanuel Tsyawo on 1/13/19.

/*
 Run linear regression for beta in Y = X*coefs + e.
 Solvers of system of normal equations available
 are 1. conjugate gradient, 2. gauss-seidel, and 3. SOR
 */
/*
 Compile using
 gcc -c linreg_test.c linreg.c matply.c solve.c utils.c read_txt.c
 gcc -o execLinreg_Test linreg_test.o linreg.o matply.o solve.o utils.o read_txt.o
 ./execLinreg_Test
 */


#include "linreg.h"
#include "utils.h"

// main function for code execution
int main( )
{
    int nrX = 1000;
    int ncX = 6;
    double *coefs, *X, *Y;
    char *datname;
    
    // allocate memory
    coefs=allocvector(ncX); Y = allocvector(nrX);
    
    datname= "dat_lreg.txt"; // name of data set in folder
    // read in data
    X=read_txt(datname, &nrX, &ncX);
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nrX, &ncX);
   
    // true parameter values to be solved for
    double trcoefs[6]={1.2,0.5,1.0,0.0,1.5,-0.5};
    
    printf("True parameter values to be solved for");
    printv(trcoefs,&ncX);
    
    printf("Calling linear regression with the conjugate gradient solver linreg_cg( ) ... \n");
    linreg_cg(Y, X, coefs, &nrX, &ncX);
    printv(coefs,&ncX);
    //printf("The solution is beta = [%.2f,%.2f]\n\n",coefs[0],coefs[1]);
    printf("Calling linear regression with the Gauss-Seidel solver linreg_gs( ) ... \n");
    linreg_gs(Y, X, coefs, &nrX, &ncX);
    printv(coefs,&ncX);
    //printf("The solution is beta = [%.2f,%.2f]\n\n",coefs[0],coefs[1]);
    printf("Calling linear regression with the SOR solver linreg_sor( ) ... \n");
    linreg_sor(Y, X, coefs, &nrX, &ncX);
    printv(coefs,&ncX);
    //printf("The solution is beta = [%.2f,%.2f]\n",coefs[0],coefs[1]);
    puts(" ");
}



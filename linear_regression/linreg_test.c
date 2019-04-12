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
#include "solve.h"

// main function for code execution
int main( )
{
    int nrX = 1000, ncX = 6, ny=1,i,j;
    double *coefs, *X, *Y, *R, z=0.0;
    char *datname;
    
    // allocate memory
    coefs=allocvector(ncX); Y = allocvector(nrX); R = allocvector(ncX*ncX);
    
    datname= "dat_lreg.txt"; // name of data set in folder
    // read in data
    X=read_txt(datname, &nrX, &ncX);
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nrX, &ncX);
   
    // true parameter values to be solved for
    double trcoefs[6]={1.2,0.5,1.0,0.0,1.5,-0.5};
    
    printf("True parameter values to be solved for");
    printm(trcoefs,&ny,&ncX);
    
    printf("Calling linear regression with the conjugate gradient solver linreg_cg( ) ... \n");
    linreg_cg(Y, X, coefs, &nrX, &ncX);
    printm(coefs,&ny,&ncX);
    
    printf("Calling linear regression with the Gauss-Seidel solver linreg_gs( ) ... \n");
    linreg_gs(Y, X, coefs, &nrX, &ncX);
    printm(coefs,&ny,&ncX);
    
    printf("Calling linear regression with the SOR solver linreg_sor( ) ... \n");
    linreg_sor(Y, X, coefs, &nrX, &ncX);
    printm(coefs,&ny,&ncX);
    
    printf("Calling linear regression via QR least squares linreg_qrc( ) ... \n");
    linreg_qrc(Y, X, coefs, &nrX, &ncX);
    printm(coefs,&ny,&ncX);
    puts(" ");
    
    
    printf("Calling linear regression via coordinate descent algorithm linreg_cord( ) ... \n");
    init_vec (coefs, &ncX, &z);
    linreg_cord(Y, X, coefs, &nrX, &ncX);
    printm(coefs,&ny,&ncX);
    puts(" ");
    
    printf("Linear regression with QR using Modified Gram-Schmidt Algorithm. \n First with version that does not overwrite X with Q\n");
    linreg_qrMGS(Y, X, coefs, R, &nrX, &ncX);
    printf("The solution is \n");
    printm(coefs,&ny,&ncX);
    
    // read in data
    X=read_txt(datname, &nrX, &ncX);
    
    // prepare data for regression
    dat_read_prep(datname, X, Y, &nrX, &ncX);

    printf("Linear regression with QR using Modified Gram-Schmidt Algorithm. \n Then with version that overwrites X with Q\n");
    linreg_qrMGS2(Y, X, coefs, R, &nrX, &ncX);
    printf("The solution is \n");
    printm(coefs,&ny,&ncX);
    
}



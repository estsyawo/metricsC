/*
  Emmanuel S. Tsyawo 
  estsyawo@temple.edu,  estsyawo@gmail.com
  February 08, 2019
  Testing functions in mat_inv.c for inverting matrices.
*/
/* Compile using
 gcc -c mat_inv_test.c mat_inv.c solve.c conjgrad.c gauss_seidel.c SOR.c matply.c
 gcc -o execMat_Inv_Test mat_inv_test.o mat_inv.o solve.o conjgrad.o gauss_seidel.o SOR.o matply.o
 ./execMat_Inv_Test
 */

#include "matply.h"

// main program
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6;
    double *Ainv;
    char solver;

    Ainv = malloc((nra*nra)*sizeof(double)); // vector to store the solution

    printf("\n");
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
        3.446,88.798,5.796,-2.039,-18.366,-15.693,
        2.686,5.796,85.911,7.203,-16.848,-16.081,
        6.798,-2.039,7.203,96.208,-1.853,-1.8,
        2.59,-18.366,-16.848,-1.853,90.32,7.953,
        -19.055,-15.693,-16.081,-1.8,7.953,109.126};

    printf("The elements of matrix A are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", A[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");

// Conjugate gradient algorithm
    printf("\n");
    printf("Calling the Conjugate gradient algorithm to solve for Z in AZ=b \n");

    solver = 'C';
    mat_inv(A, Ainv, &n, &solver);
    printf("\n");
    
    printf("The elements of matrix Ainv are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", Ainv[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");

// Gauss-Seidel Algorithm
    printf("\n");
    printf("Calling the Gauss-Seidel algorithm to solve for Z in AZ=b \n");

    solver = 'G';
    mat_inv(A, Ainv, &n, &solver);
    printf("\n");
    
    printf("The elements of matrix Ainv are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", Ainv[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");

// SOR Algorithm
    printf("\n");
    printf("Calling the SOR algorithm to solve for Z in AZ=b \n");

    solver = 'S';
    mat_inv(A, Ainv, &n, &solver);
    printf("\n");
    
    printf("The elements of matrix Ainv are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", Ainv[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");

}

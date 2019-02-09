/*
  Emmanuel S. Tsyawo 
  estsyawo@temple.edu,  estsyawo@gmail.com
  February 08, 2019
  A wrapper function for solvers of system of linear equations
*/
/* Compile using
 gcc -c solve_test.c conjgrad.c gauss_seidel.c SOR.c matply.c
 gcc -o execSolve_Test solve_test.o conjgrad.o gauss_seidel.o SOR.o matply.o
 ./execSolve_Test
 */

#include "matply.h"
#include "solve.h"
// a wrapper function for solving a system of linear equations AX=b

// prototype function declaration

void solve(double *A, double *b, double *X, int *n, char *type);

// main program
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6;
    double *Z;
    char typeC, typeG, typeS ;

    Z = malloc(nra*sizeof(double)); // vector to store the solution

    printf("\n");
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
        3.446,88.798,5.796,-2.039,-18.366,-15.693,
        2.686,5.796,85.911,7.203,-16.848,-16.081,
        6.798,-2.039,7.203,96.208,-1.853,-1.8,
        2.59,-18.366,-16.848,-1.853,90.32,7.953,
        -19.055,-15.693,-16.081,-1.8,7.953,109.126};

    double X[6] = {-1.0,-0.5,0.0,0.5,1.0,1.5};
    double b[6] = {-111.501,-90.771,-42.952,37.773,107.917,197.644};

    printf("The elements of matrix A are: \n" );

    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf(" %.3f ", A[j*nra + i]);
        }
        puts(" ");
    }
    printf("\n");

    printf("The elements of vector b are: \n");
    for (i=0; i<nra; i++) {
        printf("b[%d] = %.3f \n",i,b[i]);
    }
    printf("\n");
    printf("The vector X of values to be solved for are: \n");
    for (i=0; i<nra; i++) {
        printf("X[%d] = %.1f \n",i,X[i]);
    }

// Conjugate gradient algorithm
    printf("\n");
    printf("Calling the Conjugate gradient algorithm to solve for Z in AZ=b \n");

    typeC = 'C';
    solve(A, b, Z, &n, &typeC);
    printf("\n");
    printf("The solution is the vector: \n");

    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,Z[i]);
    }

// Gauss-Seidel Algorithm
    printf("\n");
    printf("Calling the Gauss-Seidel algorithm to solve for Z in AZ=b \n");

    typeG = 'G';
    solve(A, b, Z, &n, &typeG);
    printf("\n");
    printf("The solution is the vector: \n");

    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,Z[i]);
    }

// SOR Algorithm
    printf("\n");
    printf("Calling the SOR algorithm to solve for Z in AZ=b \n");

    typeS = 'S';
    solve(A, b, Z, &n, &typeS);
    printf("\n");
    printf("The solution is the vector: \n");

    for(i=0; i<nra; i++){
        printf(" Z[%d] = %.3f \n",i,Z[i]);
    }

}

// function for the conjugate gradient algorithm
void solve(double *A, double *b, double *X, int *n, char *type)
{
  switch(*type){
    	case 'C':
	  conjgrad(A, b, X, n);
	  break;
	case 'G':
	  gauss_seidel(A, b, X, n);
	  break;
	case 'S':
	  SOR(A, b, X, n);
	  break;
	default :
	  printf("Invalid type of solver.\n");
  }	

}

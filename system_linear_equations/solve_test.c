/*
  Emmanuel S. Tsyawo 
  estsyawo@temple.edu,  estsyawo@gmail.com
  February 08, 2019
  A wrapper function for solvers of system of linear equations
*/
/* Compile using
 gcc -c solve_test.c solve.c matply.c utils.c
 gcc -o execSolve_Test solve_test.o solve.o matply.o utils.o
 ./execSolve_Test
 */

#include "matply.h"
#include "solve.h"
#include "utils.h"
// a wrapper function for solving a system of linear equations AX=b

// main program
int main()
{
    int i=1, n = 6;
    double *Z;
    char type;

    Z = allocvector(n); // vector to store the solution

    printf("\n");
    // matrix A
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
        3.446,88.798,5.796,-2.039,-18.366,-15.693,
        2.686,5.796,85.911,7.203,-16.848,-16.081,
        6.798,-2.039,7.203,96.208,-1.853,-1.8,
        2.59,-18.366,-16.848,-1.853,90.32,7.953,
        -19.055,-15.693,-16.081,-1.8,7.953,109.126};
    // Vectors X and b
    double X[6] = {-1.0,-0.5,0.0,0.5,1.0,1.5}; // the solution vector to be solved for
    double b[6] = {-111.501,-90.771,-42.952,37.773,107.917,197.644};

    printf("The elements of matrix A are: \n" );
    printm(A,&n,&n);
    printf("\n");

    printf("The elements of vector b are: \n");
    printm(b,&i,&n);
    printf("\n");
    
    printf("The vector X of values to be solved for are: \n");
    printm(X,&i,&n);
    printf("\n");

// Conjugate gradient algorithm
    printf("Calling the Conjugate gradient algorithm to solve for Z in AZ=b \n");
    type = 'C';
    solve(A, b, Z, &n, &type);
    printf("\n");
    printf("The solution is the vector: \n");
    printm(Z,&i,&n);
    printf("\n");

// Gauss-Seidel Algorithm
    printf("Calling the Gauss-Seidel algorithm to solve for Z in AZ=b \n");
    type = 'G';
    solve(A, b, Z, &n, &type);
    printf("\n");
    printf("The solution is the vector: \n");
    printm(Z,&i,&n);
    printf("\n");
    
// SOR Algorithm
    printf("Calling the SOR algorithm to solve for Z in AZ=b \n");
    type = 'S';
    solve(A, b, Z, &n, &type);
    printf("\n");
    printf("The solution is the vector: \n");
    printm(Z,&i,&n);
    
// QR Algorithm
    printf("Calling the QR algorithm to solve for Z in AZ=b \n");
    type = 'Q';
    solve(A, b, Z, &n, &type);
    printf("\n");
    printf("The solution is the vector: \n");
    printm(Z,&i,&n);

}

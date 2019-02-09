/*
 Emmanuel S. Tsyawo
 estsyawo@temple.edu,  estsyawo@gmail.com
 February 08, 2019
 A wrapper function for solvers of system of linear equations
 */
/* Compile using
 gcc -c SOR_test.c matply.c
 gcc -o execSOR_test SOR_test.o matply.o
 ./execSOR_test
 */

#include "solve.h"

/*  
  Input:
	A - nxn matrix
	b - nx1 vector
	X - nx1 vector of unknowns to be solved for
	type - a character for type of solver. solvers currently available are
	"conjgrad", "gauss_seidel", and "SOR"
*/

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

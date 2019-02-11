// multiply two matrices
#include "matply.h"
/* xa - vectorised matrix A,
 xb - vectorised matrix B,
 xab - product AB (vectorised)
 nra - number of rows in A
 nca - number of columns in B
 ncb - number of columns in B
 */
/*compile using:
 gcc -c matply_test.c matply.c
 gcc -o execMatplytest matply_test.o matply.o
 ./execMatplytest
 */

int main()
{
    int nra, nca,nlx, ncb,i,j;
    double *xab, *txab,*xx, *xyt, *xty;
    int nlab=6;
    nra = 2; nca = 3; ncb = nra;
    nlx = 9;

    xab = malloc(nlab*sizeof(double));
    txab = malloc(nlab*sizeof(double)); //for its transpose
    xx = malloc(nlx*sizeof(double));
    xyt = malloc((nra*nra)*sizeof(double));
    xty = malloc((nca*nca)*sizeof(double));

    double xa[6] = {-1.48,-0.93,-1.02,-0.21,1.65,0.73};
    double xb[6] = {-0.23,0.58,2.42,0.02,0.34,0.16} ;

    // Take the matrix product xab = xa*xb
    matply(xa, xb, xab, &nra, &nca, &ncb);
    printf("The matrix xab = xa*xb is \n");
    for (i=0; i<nra; i++) {
        for (j=0; j<ncb; j++) {
            printf("%.2f ",xab[j*nra + i]);
        }
        puts(" ");
    }

    // taking a product xx = xa'xa and printing out results
    printf("taking a product xa'xa and printing out results\n");
    matply_sym( xa, xx, &nra, &nca);
    for (i=0; i<nca; i++) {
        for (j=0; j<nca; j++) {
            printf("%.2f ",xx[j*nca + i]);
        }
        puts(" ");
    }
    
    printf("taking a product xa'*xb and printing out results\n");
    
    matply_xty(xa, xb, xty, &nra, &nca, &nca);
    
    for (i=0; i<nca; i++) {
        for (j=0; j<nca; j++) {
            printf("%.2f ",xty[j*nca + i]);
        }
        puts(" ");
    }
    
    printf("taking a product xa*xb' and printing out results\n");
    
    matply_xyt(xa, xb, xyt, &nra, &nca, &nra);
    for (i=0; i<nra; i++) {
        for (j=0; j<nra; j++) {
            printf("%.2f ",xyt[j*nca + i]);
        }
        puts(" ");
    }

    printf("Taking the transpose of xab\n");
    trans(xab, txab, &nra, &ncb);
    for (i=0; i<nra; i++) {
        for (j=0; j<ncb; j++) {
            printf("%.2f ",txab[j*nra + i]);
        }
        puts(" ");
    }

    printf("The dot product of xa and xb = %.2f \n", dotprod(xa, xb, &nlab));
    puts(" ");
    
    printf("The maximum element of the vector xa = %.2f\n",max(xa,&nlab));
    puts("");
    
    printf("The minimum element of the vector xb = %.2f\n",min(xb,&nlab));
    puts("");
    
}// end main

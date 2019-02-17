
#include "matply.h"
#include "utils.h"
/* xa - vectorised matrix A,
 xb - vectorised matrix B,
 xab - product AB (vectorised)
 nra - number of rows in A
 nca - number of columns in B
 ncb - number of columns in B
 */
/*compile using:
 gcc -c matply_test.c matply.c utils.c
 gcc -o execMatplytest matply_test.o matply.o utils.o
 ./execMatplytest
 */

int main()
{
    int nra, nca,nlx, ncb,i,j;
    double *xab, *txab,*xx, *xyt, *xty, a, *zz;
    int nlab=6, ns = 1;
    nra = 2; nca = 3; ncb = nra;
    nlx = 9; a = 2.1;
    
    xab = allocvector(nlab);
    txab = allocvector(nlab); //for its transpose
    xx = allocvector(nlx);
    zz = allocvector(nlab*nlab);
    xyt = allocvector(nra*nra);
    xty = allocvector(nca*nca);

    double xa[6] = {-1.48,-0.93,-1.02,-0.21,1.65,0.73};
    double xb[6] = {-0.23,0.58,2.42,0.02,0.34,0.16} ;

    // Take the matrix product xab = xa*xb
    matply(xa, xb, xab, &nra, &nca, &ncb);
    printf("The matrix xab = xa*xb is \n");
    printm(xab, &nra, &ncb);

    // taking a product xx = xa'xa and printing out results
    printf("taking a product xa'xa and printing out results\n");
    matply_sym( xa, xx, &nra, &nca);
    printm(xx, &nca, &nca);
    
    printf("taking a product xa'xa with xa as row vector and printing out results\n");
    matply_sym( xa, zz, &ns, &nlab);
    printm(zz, &nlab, &nlab);
    
    printf("taking a product xa'*xb and printing out results\n");
    
    matply_xty(xa, xb, xty, &nra, &nca, &nca);
    printm(xty, &nca, &nca);
    
    
    printf("taking a product xa*xb' and printing out results\n");
    
    matply_xyt(xa, xb, xyt, &nra, &nca, &nra);
    printm(xyt, &nra, &nra);

    printf("Taking the transpose of xab\n");
    trans(xab, txab, &nra, &ncb);
    printm(txab, &nra, &ncb);

    printf("The dot product of xa and xb = %.2f \n", dotprod(xa, xb, &nlab));
    puts(" ");
    
    printf("Taking scalar-matrix/vector product ax =a*x \n");
    matply_ax(xa, &a, txab, &nlab);
    printv(txab, &nlab);
    
    printf("Adding matrix ax to itself and passing it to itself \n");
    vecadd(txab, txab, txab, &nlab);
    printm(txab, &nra, &nca);
    
    printf("Subtracting matrix ax from itself and passing it to itself \n");
    vecsub(txab, txab, txab, &nlab);
    printm(txab, &nra, &nca);
    
    printf("The maximum element of the vector xa = %.2f\n",max(xa,&nlab));
    puts("");
    
    printf("The minimum element of the vector xb = %.2f\n",min(xb,&nlab));
    puts("");
    
    printf("The p = %d norm of vector xb equals %3.3f\n",nra,norm_lp(xb,&nlab,&nra));
    
    printf("The vector xb is \n");
    printv(xb,&nlab);
    
    printf("The infinity norm of xb is = %3.3f \n",norm_max(xb,&nlab));
}// end main

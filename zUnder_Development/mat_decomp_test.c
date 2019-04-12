//
//  matrix_decompositions_test.c
//  
//
//  Created by Emmanuel Tsyawo on 3/13/19.
//
// Compile using:
/*
 gcc -c mat_decomp_test.c mat_decomp3.c matply.c utils.c read_txt.c
 gcc -o execMD mat_decomp_test.o mat_decomp3.o matply.o utils.o read_txt.o
 ./execMD
 */
#include "utils.h"
#include "matply.h"
// Prototype declaration
void qr_hh(double *A, double *Q, double *Qk, int *m, int *n);

int main()
{
    int i,j, nr=12, nc = 11, nr1=nc-1;
    double *X, *Y, *dat, *Q, *Qk;
    char *datname = "data.txt";
    
    Y = allocvector(nr); // allocate memory to Y
    X = allocvector(nr*(nc-1));
    Q = allocvector(nr*nr);
    Qk = allocvector(nr*nr);
    // read in data
    dat=read_txt(datname, &nr, &nc );
    
    for (i=0; i<nr; i++) {
        for (j=0; j<(nc-1); j++) {
            X[j*nr+i]=dat[j*nr+i];
            Y[i]=dat[(nc-1)*nr+i];
        }
    }
    
    
    printf("Print out data set \n");
    printm(dat, &nr, &nc);
    
    printf("Print out matrix X \n");
    printm(X, &nr, &nr1);

    printf("Print out vector Y \n");
    printv(Y, &nr);
    
    
    qr_hh(X, Q, Qk, &nr, &nr1);
    
    printf("Matrix Q is \n");
    printm(Q,&nr,&nr);
    
    printf("The matrix R is \n");
    printm(X,&nr,&nr1);
    
    printf("The matrix Qk is \n");
    printm(Qk,&nr,&nr);
    
    /*
    for (i=0; i<nr; i++) {
        for (j=0; j<(nc-1); j++) {
            X[j*nr+i]=dat[j*nr+i];
     //       Y[i]=dat[(nc-1)*nr+i];
        }
    }
    
    printf("The matrix R is \n");
    matply_xty(Q,X,X,&nr,&nr,&nr1);
    printm(X,&nr,&nr1);
     */
}

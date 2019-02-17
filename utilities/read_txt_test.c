//
//  read_txt_test.c
//  
//
//  Created by Emmanuel Tsyawo on 2/9/19.
//
/*
 Compile using:
 gcc -c read_txt_test.c read_txt.c utils.c
 gcc -o execReadDat read_txt_test.o read_txt.o utils.o
 ./execReadDat
 */


#include "utils.h"

int main()
{
    int i,j, nr=12, nc = 11;
    double *X, *Y;
    char *datname = "data.txt";
    
    Y = allocvector(nr); // allocate memory to Y
    
    X=read_txt(datname, &nr, &nc );
    
    printf("Print out data set \n");
    printm(X, &nr, &nc);
    
    printf("Now assign first column to vector Y and replace with 1's in X\n");
    
    dat_read_prep(datname, X, Y, &nr, &nc);
    puts("");
    printf("Print out matrix X");
    printm(X, &nr, &nc);
    
    printf("Print out vector Y \n");
    printv(Y, &nr);
}

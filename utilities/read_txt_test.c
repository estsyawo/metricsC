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
    double *data;
    char *datname = "data.txt";
    
    data=read_txt(datname, &nr, &nc );
    
    for (i=0; i<nr; i++) {
        for (j=0; j<nc; j++) {
            printf(" %.1f ",data[j*nr + i]);
        }
        puts("");
    }
}

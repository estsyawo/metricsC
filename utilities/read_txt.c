//
//  read_txt.c
//
//  Created by Emmanuel Tsyawo on 2/9/19.
//  Read a .txt data file as a matrix (nr x nc) and store in long vector dat (length nrxnc)

#include "utils.h"
/* Arguments:
    datname as "filename.txt" as a pointer
    nr - number of rows in the data file
    nc - number of columns in the data file
 Output:
    a long vector of length (nrxnc) which stores the data column major
 */

double *read_txt(char *datname, int *nr, int *nc )
{

int i, j, readerror=0, success;
double *dat;
char s[25];
FILE *datfile;
datfile=fopen(datname, "r");

    /* Allocate long vector */
dat=allocvector((*nr)*(*nc));


for (i=0;i<*nr;i++) {
    if ( readerror ) break;
    for (j=0;j<*nc;j++) {
        if (fscanf(datfile, "%s", s)==EOF) {
            readerror=1;
            printf("Not enough data\n");
            printf("i=%u j=%u\n", i, j);
            break;
        }
        dat[j*(*nr) + i] = strtod(s,NULL);
        
    }
}
    return dat;
fclose(datfile);
}

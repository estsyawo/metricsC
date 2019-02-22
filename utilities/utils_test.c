//
//  utils_test.c
//  
//
//  Created by Emmanuel Tsyawo on 2/11/19.
//

#include "utils.h"

/* Compile using:
 gcc -c utils_test.c utils.c
 gcc -o execUtils_Test utils_test.o utils.o
 ./execUtils_Test
 */

// test functions in utils.c
int main()
{
    int nra = 2, nca = 3, nlx = 6, ny=1;
    double xa[6] = {-1.48,-0.93,-1.02,-0.21,1.65,0.73};
    double xd[6] = {-1.48,0.73,-0.93,-0.21,1.65,-1.02};
    int xb[6] = {-1,2,3,26,0,8};
    double x = 0.5, mu = 0, sig = 1.0;
    
    // testing the printv() function
    printf("testing the printv() function...\n");
    printv(xa, &nlx);
    
    printf("testing the printm_int() function...\n");
    printm_int(xb, &nra, &nca);
    
    printf("testing the printm() function...\n");
    printm(xa, &nra, &nca);
    
    printf("testing the function init_vec\n initialising xa to zeros\n");
    init_vec (xa, &nlx, &mu);
    printv(xa, &nlx);
    
    printf("d_logit(0.5) = %.2f\n", d_logit(x));
    printf("p_logit(0.5) = %.2f\n", p_logit(x));
    
    printf("d_norm(0.5) = %.2f\n", d_norm(x,mu,sig));
    printf("p_norm(0.5) = %.2f\n", p_norm(x,mu,sig));
    
    printf("The vector xd of doubles before sorting \n");
    printv(xd,&nlx);
    sort(xd, nlx);
    printf("The vector xd after sorting \n");
    printv(xd,&nlx);
    
    printf("The vector xb of integers before sorting \n");
    printm_int(xb,&ny,&nlx);
    sort_int(xb, nlx);
    printf("The vector xb after sorting \n");
    printm_int(xb,&ny,&nlx);
    

}

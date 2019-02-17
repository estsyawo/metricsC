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
    int nra = 2, nca = 3, nlx = 6;
    double xa[6] = {-1.48,-0.93,-1.02,-0.21,1.65,0.73};
    double x = 0.5, mu = 0, sig = 1.0;
    
    // testing the printv() function
    printf("testing the printv() function...\n");
    printv(xa, &nlx);
    
    printf("testing the printm() function...\n");
    printm(xa, &nra, &nca);
    
    printf("testing the function init_vec\n initialising xa to zeros\n");
    init_vec (xa, &nlx, &mu);
    printv(xa, &nlx);
    
    printf("d_logit(0.5) = %.2f\n", d_logit(x));
    printf("p_logit(0.5) = %.2f\n", p_logit(x));
    
    printf("d_norm(0.5) = %.2f\n", d_norm(x,mu,sig));
    printf("p_norm(0.5) = %.2f\n", p_norm(x,mu,sig));

}

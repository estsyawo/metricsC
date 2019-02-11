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
    double x = 0.5, mu = 0, sig = 1.0;
    printf("d_logit(0.5) = %.2f\n", d_logit(x));
    printf("p_logit(0.5) = %.2f\n", p_logit(x));
    
    printf("d_norm(0.5) = %.2f\n", d_norm(x,mu,sig));
    printf("p_norm(0.5) = %.2f\n", p_norm(x,mu,sig));

}

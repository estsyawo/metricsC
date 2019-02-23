/* Compile using
 gcc -c qr_householder_test.c matply.c utils.c
 gcc -o execQR_HH qr_householder_test.o matply.o utils.o
 ./execQR_HH
 */

// QR Householder algorithm for solving the least squares problem AX=b

#include "matply.h"
#include "utils.h"
// prototype function declaration

void lsqqr( double *a, double *b, double *coef, int *nrow,  int *ncol);
// main program
int main()
{
    int i,j;
    int const nra = 6;
    int n = 6, ny = 1;
    double *Z;
    
    Z = malloc(nra*sizeof(double)); // vector to store the solution
    
    printf("\n");
    double A[36] = {87.184,3.446,2.686,6.798,2.59,-19.055,
        3.446,88.798,5.796,-2.039,-18.366,-15.693,
        2.686,5.796,85.911,7.203,-16.848,-16.081,
        6.798,-2.039,7.203,96.208,-1.853,-1.8,
        2.59,-18.366,-16.848,-1.853,90.32,7.953,
        -19.055,-15.693,-16.081,-1.8,7.953,109.126};
    
    double X[6] = {-1.0,-0.5,0.0,0.5,1.0,1.5};
    double b[6] = {-111.501,-90.771,-42.952,37.773,107.917,197.644};
    
    printf("The elements of matrix A are: \n" );
    printm(A,&n,&n);
    
    printf("\n");
    
    printf("The elements of vector b are: \n");
    printv(b,&n);
    printf("\n");
    
    printf("The vector X of values to be solved for are: \n");
    printm(X,&ny,&n);
    printf("\n");
    
    
    printf("\n");
    printf("Calling the QR Householder algorithm to solve for Z in AZ=b \n");
    
    lsqqr( A, b, Z, &n, &n);
    printf("\n");
    printf("The solution is the vector: \n");
    
    printm(Z,&ny,&n);
    printf("\n");
    
    
}

/* QR algorithm accessed from:
 Source: http://www.math.umd.edu/~mariakc/teaching-2/householder.c
 on February 22, 2019
 */

void lsqqr( double *a, double *b, double *coef, int *nrow,  int *ncol) {
    int i,m,l;
    double normx, aaux, baux;
    double aux[*ncol], u[*nrow];
    double c[*nrow]; /* weight matrix */
    
    /* put larger weights at the ends interval to reduce the error at the ends of the interval */
    /* Ax = b   <==> CAx = Cb where C = diag{c[0],c[1],...,c[nrow]}  */
    for( m=0; m<*nrow; m++ ) c[m]=1.0;
    for( m=0; m<min_2((*nrow)/2,10); m++ ) c[m]=(10-m);
    for( m=(*nrow)-1; m>max_2((*nrow)/2,(*nrow)-9); m-- ) c[m]=(9-(*nrow)+m);
    /* compute CA and Cb */
    for( m=0; m<*nrow; m++ ) {
        for( i=0; i<*ncol; i++ ) {
            (*(a+i+m*(*ncol)))*=c[m];   /* a_{m,i} = a_{m,i} * c_m */
        }
        b[m]*=c[m];   /* b_m = b_m * c_m */
    }
    
    /* Start the QR algorithm via Householder reflections */
    for( i=0; i<*ncol; i++ ) {
        /*form vector u=House(a(i:nrow,i))*/
        normx=0.0;
        for( m=0; m<*nrow-i; m++ ) {
            aaux=*(a+(i+m)*(*ncol)+i);  /* aaux = a_{m+i,i} */
            normx+=aaux*aaux;  /* normx = normx + aaux*aaux */
        }
        aaux=*(a+i*(*ncol)+i);    /* aaux = a_{i,i} */
        u[0]=aaux+sgn(aaux)*sqrt(normx);
        for( m=1; m<*nrow-i; m++ ) u[m]=*(a+(i+m)*(*ncol)+i);  /* u[m] = a_{i+m,i} */
        normx=0.0;
        for( m=0; m<*nrow-i; m++ ) normx+=u[m]*u[m];   /* normx = normx + u[m]*u[m] */
        /* compute (I-2uu^t)a(i:nrow,i:ncol) */
        for( l=0; l<*ncol-i; l++ ) {
            aux[l]=0.0;
            for( m=0; m<*nrow-i; m++ )  {
                aux[l]+=u[m]*(*(a+(i+m)*(*ncol)+i+l));  /* aux[l] = aux[l] + u[m]*a_{i+m,i+l} */
            }
            aux[l]*=(-2.0/normx);     /* aux[l] = -2*aux[l]/normx */
            for( m=0; m<*nrow-i; m++ )  {
                *(a+(i+m)*(*ncol)+i+l)+=u[m]*aux[l];      /* a_{i+m,i+l} = a_{i+m,i+l} + u[m]*aux[l] */
            }
        }
        /* compute (Q^T)*b */
        baux=0.0;
        for( m=0; m<*nrow-i; m++ ) baux+=u[m]*b[i+m];  /* baux = baux + u[m]*b[i+m] */
        
        baux*=(-2.0/normx);
        for( m=0; m<*nrow-i; m++ ) b[i+m]+=baux*u[m];  /* b[i+m] = b[i+m] + baux*u[m] */
    }
    /* find the coefficients*/
    for( i=(*ncol)-1; i>=0; i-- ) {
        baux=0.0;
        for( m=i+1; m<*ncol; m++ ) {
            baux+=(*(a+i*(*ncol)+m))*coef[m];    /* baux = baux + a_{j,m}*coef[m] */
        }
        coef[i]=(b[i]-baux)/(*(a+i*(*ncol)+i));
    }
}

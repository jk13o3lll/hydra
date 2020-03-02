#include "hydra.h"

int main(int argc, char *argv[]){
    int n = argc == 2? atoi(argv[1]) : 5;
    double *A, *x, *b, *tmp;
    // allocate
    A = (double*) malloc(n * n * sizeof(double));
    x = (double*) malloc(n * sizeof(double));
    b = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // init
    srand(time(NULL));
    // randA(n, A);
    rand_x(n * n, A, 2.0, -1.0);
    rand_x(n, b, 2.0, -1.0);
    // solve Ax = b
    gaussian(n, A, b, x);
    // bicg(n, A, b, x);
    x_minus_Ay(n, b, A, x, tmp);
    // print
    printf("A = \n"); print_A(n, A); putchar('\n');
    printf("b = \n"); print_x(n, b); putchar('\n');
    printf("x = \n"); print_x(n, x); putchar('\n');
    printf("b-Ax = \n"); print_x(n, tmp); putchar('\n');
    // release
    free(A); free(x); free(b); free(tmp);
    return 0;
}
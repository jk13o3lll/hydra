#include "hydra.h"

int main(int argc, char *argv[]){
    int n = argc == 2? atoi(argv[1]) : 3;
    double *A, *x, *b;

    // allocate
    A = (double*)malloc(n * n * sizeof(double));
    x = (double*)malloc(n * sizeof(double));
    b = (double*)malloc(n * sizeof(double));

    // initialize
    rands(n * n, A, 2.0, -1.0);
    rands(n, b, 2.0, -1.0);
    print_A(n, A, "%.8lf, ", "A = \n");
    print_x(n, b, "%.8lf, ", "b = \n");

    // // solve by bicg
    bicg(n, A, x, b);
    print_x(n, x, "%.8lf, ", "x = \n");

    // solve by gaussian
    gaussian(n, A, x, b);
    print_x(n, x, "%.8lf, ", "x = \n");

    // release
    free(A); free(x); free(b);

    return 0;
}
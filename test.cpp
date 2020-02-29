// g++ test.cpp -o test.exe -fopenmp -lpthread

#include "solver.h"

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
    randA(n, A, 2.0, -1.0);
    randx(n, b, 2.0, -1.0);
    // solve Ax = b
    // bicgstab(n, A, b, x, false, 1e-3); // sometimes cannot converge
    // bicgstab(n, A, b, x, true, 1e-3);  // preconditioned not work?
    gaussian(n, A, b, x);           // this work well
    xminusAy(n, b, A, x, tmp);
    // print
    printA("A = \n", n, A);
    printx("b = \n", n, b);
    printx("x = \n", n, x);
    printx("b-Ax = \n", n, tmp);
    // release
    free(A); free(x); free(b); free(tmp);
    return 0;
}

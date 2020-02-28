#include "bicgstab.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

int main(int argc, char *argv[]){
    // TODO: Test unpreconditioned BiCGSTAB
    // https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
    int n = argc == 2? atoi(argv[1]) : 10;
    double *A, *x, *b, *tmp;
    // allocate
    A = (double*) malloc(n * n * sizeof(double));
    x = (double*) malloc(n * sizeof(double));
    b = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // init
    srand(time(NULL));
    randA(n, A);
    // randAsparse(n, A);
    // randAsparse(n, A, n / 3, 0.7);
    randx(n, b);
    // A. bicgstab to solve
    // bicgstab(n, A, b, x, 1e-3); // return if maxIter or b-Ax approx 0
    // B. gaussian ellimination to solve
    gaussian(n, A, b, x);
    xminusAy(n, b, A, x, tmp);
    // print
    printA("A = ", n, A);
    printx("b = ", n, b);
    printx("x = ", n, x);
    printx("b - Ax = ", n, tmp);
    // release
    free(A); free(x); free(b); free(tmp);
    return 0;
}

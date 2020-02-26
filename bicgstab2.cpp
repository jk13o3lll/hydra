#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// x, y, z: vector
// A, B: matrix
// s: scalar

void zerox(int n, double *x){
    int i;
    for(i = 0; i < n; ++i)
        x[i] = 0.0;
}
void randx(int n, double *x){
    int i;
    for(i = 0; i < n; ++i)
        x[i] = (double) rand() / RAND_MAX;
}
void randA(int n, double *A){
    int i, j;
    for(i = 0; i < n; ++i)
        for(j = 0; j < n; ++j)
            A[i*n+j] = (double) rand() / RAND_MAX;
}
void xTy(int n, double *x, double *y, double &res){
    int i;
    res = 0.0;
    for(i = 0; i < n; ++i)
        res += x[i] * y[i];
}
void Ax(int n, double *A, double *x, double *res){
    int i, j;
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[i*n+j] * x[j];
    }
}
void ATx(int n, double *A, double *x, double *res){
    int i, j;
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[j*n+i] * x[j];
    }
}
void xplusy(int n, double *x, double *y, double *res){
    int i;
    for(i = 0; i < n; ++i)
        res[i] = x[i] + y[i];
}
void xminusy(int n, double *x, double *y, double *res){
    int i;
    for(i = 0; i < n; ++i)
        res[i] = x[i] - y[i];
}
void xplussy(int n, double *x, double s, double *y, double *res){
    int i;
    for(i = 0; i < n; ++i)
        res[i] = x[i] + s * y[i];
}
void xminussy(int n, double *x, double s, double *y, double *res){
    int i;
    for(i = 0; i < n; ++i)
        res[i] = x[i] - s * y[i];
}

int main(int argc, char *argv[]){
    // TODO: Try unpreconditioned BiCGSTAB
    // https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method

    return 0;
}

// int main(int argc, char *argv[]){
//     double *x, *y, *z;
//     double *A;
//     double tmp;
//     int n = 3;
//     int i, j;
    
//     // allocate
//     x = (double*) malloc(n * sizeof(double));
//     y = (double*) malloc(n * sizeof(double));
//     z = (double*) malloc(n * sizeof(double));
//     A = (double*) malloc(n * n * sizeof(double));

//     // init
//     srand(time(NULL));
//     randx(n, x);
//     randx(n, y);
//     randA(n, A);

//     // dot
//     xTy(n, x, y, tmp);
//     printf("x: "); for(i = 0; i < n; ++i) printf("%.2lf, ", x[i]); putchar('\n');
//     printf("y: "); for(i = 0; i < n; ++i) printf("%.2lf, ", y[i]); putchar('\n');
//     printf("xTy: %.2lf\n", tmp);

//     // matrix vector
//     Ax(n, A, x, z);
//     printf("A:\n");
//     for(i = 0; i < n; ++i){ for(j = 0; j < n; ++j) printf("%.2lf, ", A[i*n+j]); putchar('\n'); }
//     printf("Ax:\n");
//     for(i = 0; i < n; ++i) printf("%.2lf, ", z[i]);
//     putchar('\n');

//     // vector addition, subtraction
//     xplusy(n, x, y, z);
//     printf("x+y: "); for(i = 0; i < n; ++i) printf("%.2lf, ", z[i]); putchar('\n');
//     xminusy(n, x, y, z);
//     printf("x+y: "); for(i = 0; i < n; ++i) printf("%.2lf, ", z[i]); putchar('\n');
//     xplussy(n, x, 0.5, y, z);
//     printf("x+y: "); for(i = 0; i < n; ++i) printf("%.2lf, ", z[i]); putchar('\n');
//     xplussy(n, x, 0.5, y, z);
//     printf("x+y: "); for(i = 0; i < n; ++i) printf("%.2lf, ", z[i]); putchar('\n');

//     // release
//     free(x);
//     free(y);
//     free(z);
//     free(A);
//     return 0;
// }
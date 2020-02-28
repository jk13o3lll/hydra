#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

// x, y, z: vector
// A, B: matrix
// s: scalar
bool allwithin(int n, double *x, double s){
    int i;
    for(i = 0; i < n; ++i)
        if(x[i] > s || x[i] < -s)
            break;
    return i == n; // all x <= abs(s)
}
void zerox(int n, double *x){
    int i;
    for(i = 0; i < n; ++i)
        x[i] = 0.0;
}
void randx(int n, double *x){ // 0 ~ 1
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
double xTy(int n, double *x, double *y){
    int i;
    double res = 0.0;
    for(i = 0; i < n; ++i)
        res += x[i] * y[i];
    return res;
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
void xminusAy(int n, double *x, double *A, double *y, double *res){
    int i, j;
    for(i = 0; i < n; ++i){
        res[i] = x[i];
        for(j = 0; j < n; ++j)
            res[i] -= A[i*n+j] * y[j];
    }
}
void printx(const char *name, int n, double *x){
    int i;
    puts(name);
    for(i = 0; i < n; ++i)
        printf("%.3lf ", x[i]);
    putchar('\n');
}
void printA(const char *name, int n, double *A){
    int i, j;
    puts(name);
    for(i = 0; i < n; ++i){
        for(j = 0; j < n; ++j)
            printf("%.2lf ", A[i*n+j]);
        putchar('\n');
    }
}

void bicgstab(int n, double *A, double *b, double *x, int maxIter = 10000, double tol = 1e-6){
    double *r0, *r, *p, *v, *s, *t, *tmp; // r0 is actuall r0 transpose
    double rhoi, rhoj, beta, alpha, w;  // j = i + 1
    int i, j;
    // allocate
    r0 = (double*) malloc(n * sizeof(double));
    r = (double*) malloc(n * sizeof(double));
    p = (double*) malloc(n * sizeof(double));
    v = (double*) malloc(n * sizeof(double));
    s = (double*) malloc(n * sizeof(double));
    t = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // init
    randx(n, x); // zerox(n, x);
    xminusAy(n, b, A, x, r0);
    memcpy(r, r0, n * sizeof(double)); // or add some rand for n-1 dim
    zerox(n, p);
    zerox(n, v);
    rhoi = alpha = w = 1.0;
    // start iteration
    for(i = 0; i < maxIter; ++i){
        // stage 1
        rhoj = xTy(n, r0, r);
        beta = (rhoj / rhoi) * (alpha / w);
        xminussy(n, p, w, v, tmp);
        xplussy(n, r, beta, tmp, p);
        Ax(n, A, p, v);
        alpha = rhoj / xTy(n, r0, v);
        xplussy(n, x, alpha, p, x);
        // if (x is accurate enough) then quit
        xminusAy(n, b, A, x, tmp);
        if(allwithin(n, tmp, tol)) break;
        // stage 2
        xminussy(n, r, alpha, v, s);
        Ax(n, A, s, t);
        w = xTy(n, t, s) / xTy(n, t, t);
        xplussy(n, x, w, s, x);
        // if (x is accurate enough) then quit
        xminusAy(n, b, A, x, tmp);
        if(allwithin(n, tmp, tol)) break;


        // update
        xminussy(n, s, w, t, r);
        rhoi = rhoj;
    }
    if(i < maxIter)
        printf("Converge at %d iteration.\n", i+1);
    else
        printf("Cannot converge.\n");
    // release
    free(r0); free(r); free(p); free(v); free(s); free(t); free(tmp);
}
int main(int argc, char *argv[]){
    // TODO: Test unpreconditioned BiCGSTAB
    // https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
    const int n = 10;
    double *A, *x, *b, *tmp;
    // allocate
    A = (double*) malloc(n * n * sizeof(double));
    x = (double*) malloc(n * sizeof(double));
    b = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // init
    srand(time(NULL));
    randA(n, A);
    randx(n, b);
    // bicgstab to solve
    bicgstab(n, A, b, x); // return if maxIter or b-Ax approx 0
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

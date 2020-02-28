#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

// x, y, z: vector; A, B: matrix; s: scalar
int allwithin(int n, double *x, double s){
    int i;
    for(i = 0; i < n; ++i)
        if(isnan(x[i]) || x[i] > s || x[i] < -s)
            break;
    if(i == n)              return 1; // all within
    else if(isnan(x[i]))    return 2; // some are nan
    else                    return 0; // some are not within
}
void zerox(int n, double *x){
    int i;
    for(i = 0; i < n; ++i)
        x[i] = 0.0;
}
void randx(int n, double *x, double a = 1.0, double b = 0.0){ // a*x+b (x: 0 ~ 1)
    int i;
    for(i = 0; i < n; ++i)
        x[i] = (double) rand() / RAND_MAX * a + b;
}
void randA(int n, double *A){
    int i, j;
    for(i = 0; i < n; ++i)
        for(j = 0; j < n; ++j)
            A[i*n+j] = (double) rand() / RAND_MAX;
}
void randAsparse(int n, double *A, int a = 0, double p = 0.9, bool permute = true){
    // a: cover at most how many neighbor variable (backward) (0<a<n)
    // p: probability to set random value to that neighbor (0<p<1)
    int i, j, ii;
    double tmp;
    // adjust parameters
    if(a == 0)
        a = n > 12? n/4 : 3;
    // at least non zero random on diagonal
    for(i = 0; i < n; ++i)
        for(j = 0; j < n; ++j)
            A[i*n+j] = i!=j? 0.0 : ((double)rand()/RAND_MAX*2.0-1.0);
    // distribute some random values
    for(i = 0; i < n; ++i)
        for(j = 1; j < a; ++j)
            if((double)rand() < p * RAND_MAX)
                A[i*n + (i+j)%n] = (double)rand()/RAND_MAX*2.0-1.0;
    // some column and row permutation (change order of variables and order of equations)
    if(permute)
        for(i = 0; i < n; ++i){
            // row permutation
            ii = rand() % n;
            for(j = 0; j < n; ++j){ tmp = A[i*n+j]; A[i*n+j] = A[ii*n+j]; A[ii*n+j] = tmp; }
            // column permutation
            ii = rand() % n;
            for(j = 0; j < n; ++j){ tmp = A[j*n+i]; A[j*n+i] = A[j*n+ii]; A[j*n+ii] = tmp; }
        }
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
double sx(int n, double s, double *x, double *res){
    int i;
    for(i = 0; i < n; ++i)
        res[i] = s * x[i];
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

// https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
// good for large and sparse, may not converge for large and dense systems
void bicgstab(int n, double *A, double *b, double *x, int maxIter = 1000, double tol = 1e-6, int maxAttemp = 1000){
    double *r0, *r, *p, *v, *s, *t, *tmp; // r0 is actuall r0 transpose
    double rhoi, rhoj, beta, alpha, w;  // j = i + 1
    int i, j, res;

    // adjust max iteration based on matrix size
    if(10 * n > maxIter)
        maxIter = 10 * n;
    if(100 * n > maxAttemp)
        maxAttemp = 100 * n;
    // allocate
    r0 = (double*) malloc(n * sizeof(double));
    r = (double*) malloc(n * sizeof(double));
    p = (double*) malloc(n * sizeof(double));
    v = (double*) malloc(n * sizeof(double));
    s = (double*) malloc(n * sizeof(double));
    t = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // attempts
    for(i = 0; i < maxAttemp; ++i){
        printf("Attempt %d:\n", i + 1);
        // init
        if(i == 0)                      zerox(n, x);
        else if(i * 3 < maxAttemp)      randx(n, x, 2.0, -1.0);
        else if(i * 3 < maxAttemp * 2)  randx(n, x, 200.0, -100.0);
        else                            randx(n, x, 0.02, -0.01);
        xminusAy(n, b, A, x, r0);
        memcpy(r, r0, n * sizeof(double)); // or add some rand for n-1 dim
        zerox(n, p);
        zerox(n, v);
        rhoi = alpha = w = 1.0;
        // start iteration
        for(j = 0; j < maxIter; ++j){
            // stage 1
            rhoj = xTy(n, r0, r);
            if((rhoj < 1e-12 && rhoj > -1e-12) || (rhoj < -1e12 || rhoj > 1e12))
                rhoj = rhoi; // avoid rho become too small or too large
            beta = (rhoj / rhoi) * (alpha / w);
            xminussy(n, p, w, v, tmp);
            xplussy(n, r, beta, tmp, p);
            Ax(n, A, p, v);
            alpha = rhoj / xTy(n, r0, v);
            xplussy(n, x, alpha, p, x);
            // if x is accurate enough, then quit
            xminusAy(n, b, A, x, tmp);
            if((res = allwithin(n, tmp, tol)) > 0) break;
            // stage 2
            xminussy(n, r, alpha, v, s);
            Ax(n, A, s, t);
            w = xTy(n, t, s) / xTy(n, t, t);
            xplussy(n, x, w, s, x);
            // if x is accurate enough, then quit
            xminusAy(n, b, A, x, tmp);
            if((res = allwithin(n, tmp, tol)) > 0) break;
            // update
            xminussy(n, s, w, t, r);
            rhoi = rhoj;
        }
        if(j == maxIter)
            printf("Cannot converge.\n");
        else if(res == 1){
            printf("Converge at %d iteration.\n", i+1);
            break;
        }
        else if(res == 2)
            printf("Detect NaN.\n");
    }
    // release
    free(r0); free(r); free(p); free(v); free(s); free(t); free(tmp);
}


// for large and dense, try Gauss-Jordan?
// https://en.wikipedia.org/wiki/Gaussian_elimination
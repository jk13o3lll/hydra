// TODO: still has bug, check fundamental fnctions for vector and matrix calculation is correct first

#include <stdio.h>
#include <stdlib.h> // malloc, free, srand, rand
#include <string.h> // memcopy
#include <time.h>   // time

// s for scalar;  x, y for vector;  A for matrix

// Use memcpy to copy
// memcpy(y, x, n * sizeof(double)); // y = x

// return scalar
double xTy(double *x, double *y, int n){
    double res = 0.0;
    int i;
    for(i = 0; i < n; ++i)
        res += x[i] * y[i];
    return res;
}
double xTAy(double *x, double *A, double *y, int n){
    double res = 0.0;
    int i, j;
    for(i = 0; i < n; ++i)
        for(j = 0; j < n; ++j)
            res += x[i] * A[i*n+j] * y[j];
    return res;
}

// if res != NULL, use memory already allocated
double* sx(double s, double *x, int n, double *res = NULL){
    int i;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i)
        res[i] = s * x[i];
    return res;
}
double* xAddy(double *x, double *y, int n, double *res = NULL){ // res = x for inplace
    int i;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i)
        res[i] = x[i] + y[i];
    return res;
}
double* xSuby(double *x, double *y, int n, double *res = NULL){
    int i;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i)
        res[i] = x[i] - y[i];
    return res;
}
double* xAddsy(double *x, double s, double *y, int n, double *res = NULL){ // res = x for inplace
    int i;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i)
        res[i] = x[i] + s * y[i];
    return res;
}
double* xSubsy(double *x, double s, double *y, int n, double *res = NULL){
    int i;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i)
        res[i] = x[i] - s * y[i];
    return res;
}
double* Ax(double *A, double *x, int n, double *res = NULL){
    int i, j;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[i*n+j] * x[j];
    }
    return res;
}
double* ATx(double *A, double *x, int n, double *res = NULL){
    int i, j;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[j*n+i] * x[j];
    }
    return res;
}
double* sAx(double s, double *A, double *x, int n, double *res = NULL){
    int i, j;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[i*n+j] * x[j];
        res[i] *= s;
    }
    return res;
}
double* sATx(double s, double *A, double *x, int n, double *res = NULL){
    int i, j;
    if(res == NULL)
        res = (double*)malloc(n * sizeof(double));
    for(i = 0; i < n; ++i){
        res[i] = 0.0;
        for(j = 0; j < n; ++j)
            res[i] += A[j*n+i] * x[j];
        res[i] *= s;
    }
    return res;
}

// https://docs.nvidia.com/cuda/incomplete-lu-cholesky/index.html
// https://blog.csdn.net/langb2014/article/details/51348673
// https://en.wikipedia.org/wiki/Biconjugate_gradient_stabilized_method
// https://github.com/scipy/scipy/blob/v0.14.0/scipy/sparse/linalg/isolve/iterative.py
// https://github.com/dimikout3/Parallel-BiCGSTAB-/blob/master/serial.cpp
double* bicgstab(double *A, double *b, int n, double *x = NULL){
    double *r0_, *r, *v, *p, *s, *t, *tmp; // j = i - 1, x_ = xT
    double rhoi, rhoj, alpha, beta, w, err;
    double tol = 1e-6;
    int i;

    // allocate
    if(x == NULL)
        x = (double*)malloc(n * sizeof(double));
    r0_ = (double*)malloc(n * sizeof(double));
    r = (double*)malloc(n * sizeof(double));
    v = (double*)malloc(n * sizeof(double));
    p = (double*)malloc(n * sizeof(double));
    s = (double*)malloc(n * sizeof(double));
    t = (double*)malloc(n * sizeof(double));
    tmp = (double*)malloc(n * sizeof(double));
    // init
    for(i = 0; i < n; ++i) x[i] = (double)rand() / RAND_MAX * 2.0 - 1.0; // -1~+1
    xSuby(b, Ax(A, x, n, tmp), n, r);
    memcpy(r0_, r, n * sizeof(double));
    rhoj = alpha = w = 1.0;
    for(i = 0; i < n; ++i) v[i] = p[i] = 0.0;
    // bicgstab
    while(1){
        rhoi = xTy(r0_, r, n);
        beta = (rhoi / rhoj) * (alpha / w);
        xAddsy(r, beta, xSubsy(p, w, v, n, tmp), n, p);
        Ax(A, p, n, v);
        alpha = rhoi / xTy(r0_, v, n);
        xAddsy(x, alpha, p, n, x);
        
        // if x is accurate enough, then quit
        xSuby(b, Ax(A, x, n, tmp), n, tmp);
        err = xTy(tmp, tmp, n);
        printf("%lf ", err);
        if(err < tol) break;

        xSubsy(r, alpha, v, n, s);
        Ax(A, s, n, t);
        w = xTy(t, s, n) / xTy(t, t, n);
        xAddsy(x, w, s, n, x);

        // if x is accurate enough, then quit
        xSuby(b, Ax(A, x, n, tmp), n, tmp);
        err = xTy(tmp, tmp, n);
        printf("%lf\n", err);
        if(err < tol) break;

        xSubsy(s, w, t, 3, r);
        rhoj = rhoi;
        system("pause");
    }

    // release
    free(r0_); free(r); free(v); free(p); free(s); free(t); free(tmp);
    return x;
}



int main(int argc, char *argv[]){
    double A[16] = {
        0.4, 1.1, 0.0, 0.0,
        0.4, 1.8, 0.0, 0.1,
        0.0, 3.0, -2.7, 0.0,
        0.13, 0.0, -3.2, 1.6
    };    
    double b[4] = {
        0.0, 1.0, 2.7, -0.3
    };
    double x[4];
    int i;

    srand(time(NULL));
    bicgstab(A, b, 4, x);

    for(i = 0; i < 4; ++i)
        printf("%lf ", x[i]);
    putchar('\n');

    return 0;
}

/*
import numpy as np

A = np.array([[0.4, 1.1, 0, 0], [0.4, 1.8, 0, 0.1], [0, 3, -2.7, 0], [0.13, 0, -3.2, 1.6]])
b = np.array([[0], [1], [2.7], [-0.3]])

x = np.matmul(np.linalg.inv(A), b)

>>> x
array([[-3.54825666],
       [ 1.29027515],
       [ 0.43363905],
       [ 0.96807396]])
*/
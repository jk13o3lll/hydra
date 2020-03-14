#ifndef __HYDRA_H__
#define __HYDRA_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <queue>

#define ISZERO(x) (fabs((x)) < 1e-12)

// n: size; s: scale;  x, y, z: vector; A: matrix (square); T: transpose;
int all_within(int n, double *x, double lowerb, double upperb);
void rands(int n, double *x, double a = 1.f, double b = 0.f); // x = [(0~1)*a+b, ...]
void zeros(int n, double *x); // x = [0, 0, ...]
void ones(int n, double *x); // x = [1, 1, ...]
void x_plus_sy(int n, double *x, double s, double *y, double *z); // z = x + sy
void x_minus_Ay(int n, double *x, double *A, double *y, double *z); // z = x - Ay
void x_minus_ATy(int n, double *x, double *A, double *y, double *z); // z = x - ATy
void x_minus_sAy(int n, double *x, double s, double *A, double *y, double *z); // z = x - sAy
void x_minus_sATy(int n, double *x, double s, double *A, double *y, double *z); // z = x - sATy
double xTy(int n, double *x, double *y); // return = xTy
double xTAy(int n, double *x, double *A, double *y); // return = xTAy
void print_x(int n, double *x, const char *fmt = "%.5lf ", const char *prefix = NULL);
void print_A(int n, double *A, const char *fmt = "%.5lf ", const char *prefix = NULL);

// biconjugate gradient (solve x, for Ax = b)
void bicg(int n, double *A, double *x, double *b,
    int maxattempts = 100, int maxiters = 1000, double tol = 1e-6);
// gaussian ellimination (LU with partial pivoting)
void gaussian(int n, double *A, double *x, double *b);

// void load_network();
// void get_equations();

// use only 1/10 original step
// with higer dim, you should not use small abs(tol)
// void solve(double step = 0.1, double tol = 1e-6);

// #ifdef USE_CUDA

// #endif // USE_CUDA

#endif
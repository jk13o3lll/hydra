#include <stdio.h>
#include <stdlib.h> // rand
#include <time.h>

// indexing
#define RC2I(r,c,cols) ((r)*(cols)+(c))
// terminate criteria
#define K 0.5 // 20% of step
#define MAX_ITER 100000
#define TOL 1e-5
// dimensions
#define N 4

bool bicg(const double *A, double *x, const double *b, const int n){
    double *xh, *r, *rh, *p, *ph; // h = hat
    double *r1, *rh1, *ptmp;
    double alpha, beta;
    double tmp;
    int i, j, k;

    xh = (double*)malloc(n * sizeof(double));
    r = (double*)malloc(n * sizeof(double));
    rh = (double*)malloc(n * sizeof(double));
    p = (double*)malloc(n * sizeof(double));
    ph = (double*)malloc(n * sizeof(double));
    r1 = (double*)malloc(n * sizeof(double));
    rh1 = (double*)malloc(n * sizeof(double));

    // initial guess
    for(i = 0; i < n; ++i){
        x[i] = (double)rand() / RAND_MAX * 2 - 1.0;
        // xh[i] = (double)rand() / RAND_MAX * 2 - 1.0;
        xh[i] = x[i];
    }
    // initialize others
    for(i = 0; i < n; ++i){
        r[i] = rh[i] = b[i];
        for(j = 0; j < n; ++j){
            r[i] -= A[RC2I(i, j, n)] * x[j];
            rh[i] -= A[RC2I(j, i, n)] * xh[j];
        }
        p[i] = r[i];
        ph[i] = rh[i];
    }
    // start iterate
    for(i = 0; i < MAX_ITER; ++i){
        // get alpha
        alpha = tmp = 0.0;
        for(j = 0; j < n; ++j){
            alpha += r[j] * rh[j];
            for(k = 0; k < n; ++k)
                tmp += ph[j] * A[RC2I(j, k, n)] * p[k];
        }
        alpha /= tmp;
        alpha *= K;
        // update x, get new r
        for(j = 0; j < n; ++j){
            x[j] += alpha * p[j];
            xh[j] += alpha * ph[j];
            r1[j] = r[j];
            rh1[j] = rh[j];
            for(k = 0; k < n; ++k){
                r1[j] -= alpha * A[RC2I(j, k, n)] * p[k];
                rh1[j] -= alpha * A[RC2I(k, j, n)] * ph[k];
            }
        }
        // // if element abs(Ax - b) all < tol, then stop
        // for(j = 0; j < n; ++j){
        //     tmp = 0.0;
        //     for(k = 0; k < n; ++k)
        //         tmp += A[RC2I(j, k, n)] * x[k];
        //     if(tmp > TOL || tmp < -TOL) break;
        // }
        // if(tmp > -TOL && tmp < TOL) break; // all entry within tol
        // if(alpha < TOL) break;

        // get beta
        beta = tmp = 0.0;
        for(j = 0; j < n; ++j){
            beta += r1[j] * rh1[j];
            tmp += r[j] * rh[j];
        }
        beta /= tmp;
        // update p
        for(j = 0; j < n; ++j){
            p[j] = r1[j] + beta * p[j];
            ph[j] = rh1[j] + beta * ph[j];
        }
        // update r by swap
        ptmp = r, r = r1, r1 = ptmp;
        ptmp = rh, rh = rh1, rh1 = ptmp;

        printf("%lf %lf\n", alpha, beta);
        // system("pause");
    }

    free(xh);
    free(r);
    free(rh);
    free(p);
    free(ph);
    free(r1);
    free(rh1);
    // return i < MAX_ITER;
    return true;
}

bool bicgstab(const double *A, double *x, const double *b, const int n){
    return true;
}

int main(int argc, char *argv[]){
    double A[N*N] = {
        0.4, 1.1, 0.0, 0.0,
        0.4, 1.8, 0.0, 0.1,
        0.0, 3.0, -2.7, 0.0,
        0.13, 0.0, -3.2, 1.6
    };    
    double b[N] = {
        0.0, 1.0, 2.7, -0.3
    };
    double x[N];
    int i;
    bool ret;

    srand(time(NULL));
    ret = bicg(A, x, b, N);

    if(ret){
        for(i = 0; i < N; ++i)
            printf("%lf ", x[i]);
        putchar('\n');
    }
    else
        puts("Cannot converge.");


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
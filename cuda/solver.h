#ifndef __SOLVER_H__
#define __SOLVER_H__

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>

// s: scalar, x, y: vector, A: matrix, T: transpose
int all_zero(int n, double *x, double tol = 1e-6){
    int i;
    for(i = 0; i < n; ++i)
        if(isnan(x[i]) || x[i] > tol || x[i] < -tol)
            break;
    if(i == n)              return 1;
    else if(isnan(x[i]))    return 2;
    else                    return 0;
}
void rand_x(int n, double *x, double a = 1.0, double b = 0.0){
    for(int i = 0; i < n; ++i)
        x[i] = (double) rand() / RAND_MAX * a + b;
}
void zeros_x(int n, double *x){
    memset(x, 0, n * sizeof(double));
}
void ones_x(int n, double *x){
    for(int i = 0; i < n; ++i)
        x[i] = 1.0;
}
void x_plus_sy(int n, double *x, double s, double *y, double *res){
    for(int i = 0; i < n; ++i)
        res[i] = x[i] + s * y[i];
}
void x_minus_sy(int n, double *x, double s, double *y, double *res){
    for(int i = 0; i < n; ++i)
        res[i] = x[i] - s * y[i];
}
void x_minus_Ay(int n, double *x, double *A, double *y, double *res){
    memcpy(res, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] -= A[i*n+j] * y[j];
}
void x_minus_ATy(int n, double *x, double *A, double *y, double *res){
    memcpy(res, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] -= A[j*n+i] * y[j];
}
void x_minus_sAy(int n, double *x, double s, double *A, double *y, double *res){
    memcpy(res, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] -= s * A[i*n+j] * y[j];
}
void x_minus_sATy(int n, double *x, double s, double *A, double *y, double *res){
    memcpy(res, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] -= s * A[j*n+i] * y[j];
}
double xTy(int n, double *x, double *y){
    double res = 0.0;
    for(int i = 0; i < n; ++i)
        res += x[i] * y[i];
    return res;
}
double xTAy(int n, double *x, double *A, double *y){
    double res = 0.0;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res += x[i] * A[i*n+j] * y[j];
    return res;
}
double Ax(int n, double *A, double *x, double *res){
    memset(res, 0, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] += A[i*n+j] * x[j];
}
void print_x(int n, double *x, const char *fmt = "%.3lf "){
    for(int i = 0; i < n; ++i)
        printf(fmt, x[i]);
}
void print_A(int n, double *A, const char *fmt = "%.3lf "){
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            printf(fmt, A[i*n+j]);
        putchar('\n');
    }
}

// bi-conjugate gradient (unstable algorithm, but good for parallelize)
// bicg has some convergence issue on large matrix?
void bicg(int n, double *A, double *b, double *x, int maxattempts = 100, int maxiter = 1000, double tol = 1e-6){
    int i, j, k;
    double *x_, *ri, *rj, *ri_, *rj_, *p, *p_, *tmp; // j = i+1
    double alpha, beta;

    // get parameters
    if(n * 10 > maxiter) maxiter = n * 10;

    // allocate
    x_ = (double*) malloc(n * sizeof(double));
    ri = (double*) malloc(n * sizeof(double));
    rj = (double*) malloc(n * sizeof(double));
    ri_ = (double*) malloc(n * sizeof(double));
    rj_ = (double*) malloc(n * sizeof(double));
    p = (double*) malloc(n * sizeof(double));
    p_ = (double*) malloc(n * sizeof(double));
    tmp = (double*) malloc(n * sizeof(double));
    // try several times (different ic)
    for(i = 0; i < maxattempts; ++i){
        printf("Attempt %d\n", i + 1);
        // init
        if(i == 0)      zeros_x(n, x);
        else if(i == 2) ones_x(n, x);
        else            rand_x(n, x, 2.0, -1.0);
        memcpy(x_, x, n * sizeof(double));
        x_minus_Ay(n, b, A, x, ri);
        x_minus_ATy(n, b, A, x, ri_);
        memcpy(p, ri, n * sizeof(double));
        memcpy(p_, ri_, n * sizeof(double));
        // start iteration
        for(j = 0; j < maxiter; ++j){
            alpha = xTy(n, ri_, ri) / xTAy(n, p_, A, p);
            x_plus_sy(n, x, alpha, p, x);
            x_plus_sy(n, x_, alpha, p_, x_);
            x_minus_sAy(n, ri, alpha, A, p, rj);
            x_minus_sATy(n, ri_, alpha, A, p_, rj_);
            beta = xTy(n, rj_, rj) / xTy(n, ri_, ri);
            x_plus_sy(n, rj, beta, p, p);
            x_plus_sy(n, rj_, beta, p_, p_);
            memcpy(ri, rj, n * sizeof(double));
            memcpy(ri_, rj_, n * sizeof(double));
            // check
            x_minus_Ay(n, b, A, x, tmp); // residual
            if((k = all_zero(n, tmp, tol)) > 0) break;
            // debug
            // printf("x = "); print_x(n, x); putchar('\n');
            // printf("b - Ax = "); print_x(n, tmp); getchar();
        }
        // check
        if(j == maxiter)
            puts("Cannot converge.");
        else if(k == 2)
            puts("Detect NaN.");
        else{
            printf("Converge at iteration %d.\n", j + 1);
            break;
        }
    }
    // release
    free(x_); free(ri); free(rj); free(ri_); free(rj_); free(p); free(p_); free(tmp);
}

// cannot converge well without preconditioned and large matrix
// void bicgstab(int n, double *A, double *b, double *x, int maxattempts = 100, int maxiter = 1000, double tol = 1e-6){
//     int i, j, k;
//     double *r0, *r, *p, *v, *s, *t, *tmp;
//     double rhoi, rhoj, alpha, beta, omega; // j = i+1
//     // get parameters
//     if(n * 10 > maxiter) maxiter = n * 10;
//     // allocate
//     r0 = (double*) malloc(n * sizeof(double));
//     r = (double*) malloc(n * sizeof(double));
//     p = (double*) malloc(n * sizeof(double));
//     v = (double*) malloc(n * sizeof(double));
//     s = (double*) malloc(n * sizeof(double));
//     t = (double*) malloc(n * sizeof(double));
//     tmp = (double*) malloc(n * sizeof(double));
//     // several attempts with different initial guess
//     for(i = 0; i < maxattempts; ++i){
//         printf("Attempt %d\n", i+1);
//         // init
//         if(i == 0)      zeros_x(n, x);
//         else if(i == 1) ones_x(n, x);
//         else            rand_x(n, x, 2.0, -1.0);
//         x_minus_Ay(n, b, A, x, r);
//         memcpy(r0, r, n * sizeof(double));
//         rhoi = alpha = omega = 1.0;
//         zeros_x(n, p);
//         zeros_x(n, v);
//         for(j = 0; j < maxiter; ++j){
//             // stage 1
//             rhoj = xTy(n, r0, r);
//             beta = (rhoj * alpha) / (rhoi * omega);
//             x_minus_sy(n, p, omega, v, tmp);
//             x_plus_sy(n, r, beta, tmp, p);
//             Ax(n, A, p, v);
//             alpha = rhoj / xTy(n, r0, v);
//             x_plus_sy(n, x, alpha, p, x);
//             // check
//             x_minus_Ay(n, b, A, x, tmp);
//             if((k = all_zero(n, tmp, tol)) > 0) break;
//             // stage2
//             x_minus_sy(n, r, alpha, v, s);
//             Ax(n, A, s, t);
//             omega = xTy(n, t, x) / xTy(n, t, t);
//             x_plus_sy(n, x, omega, s, x);
//             // check
//             x_minus_Ay(n, b, A, x, tmp);
//             if((k = all_zero(n, tmp, tol)) > 0) break;
//             // update
//             x_minus_sy(n, s, omega, t, r);
//             rhoi = rhoj;
//         }
//         // check
//         if(j == maxiter)
//             puts("Cannot converge.");
//         else if(k == 2) // break due to NaN
//             puts("Detect NaN.");
//         else{ // converge
//             printf("Converge at iteration %d\n", j + 1);
//             break;
//         }
//     }
//     // release
//     free(r0); free(r); free(p); free(v); free(s); free(t); free(tmp);
// }

void gaussian(int n, double *A, double *b, double *x){
    double *A_, *b_;
    int pr, pc, pi; // pivot row, pivot column, pivot index
    int i, j, ii, jj;
    double gmax, tmp;

    // allocate
    A_ = (double*) malloc(n * n * sizeof(double));
    b_ = (double*) malloc(n * sizeof(double));
    // init
    memcpy(A_, A, n * n * sizeof(double));
    memcpy(b_, b, n * sizeof(double)); // x as temporary b
    // To row encholon form
    for(pr = pc = 0; pr < n && pc < n; ++pc){
        // find max abs
        pi = pr, gmax = fabs(A_[pr*n+pc]);
        for(i = pr+1; i < n; ++i)
            if((tmp = fabs(A_[i*n+pc])) > gmax)
                pi = i, gmax = tmp;
        // start ellimination
        if(gmax < 1e-12)        // pivot is zeri
            A_[pr*n+pc] = NAN;  // for convenience later
        else{
            if(pi != pr){ // swap row
                memcpy(x, A_+pr*n, n * sizeof(double));
                memcpy(A_+pr*n, A_+pi*n, n * sizeof(double));
                memcpy(A_+pi*n, x, n * sizeof(double));
                tmp = b_[pr], b_[pr] = b_[pi], b_[pi] = tmp;
            }
            for(i = pr + 1; i < n; ++i){ // elliminate
                tmp = -A_[i*n+pc] / A_[pr*n+pc];
                A_[i*n+pc] = 0.0;
                #pragma omp parallel for
                for(j = pc + 1; j < n; ++j)
                    A_[i*n+j] += A_[pr*n+j] * tmp;
                b_[i] += b_[pr] * tmp;
            }
            ++pr;
        }
    }
    // retrieve x
    memcpy(x, b_, n * sizeof(double));
    for(pr = pc = n - 1; pr >= 0 && pc >=0; ){
        if(pr + 1 < n && isnan(A_[(pr+1)*n+pc])){
            printf("Has multiple solutions for x[%d]\n", pc);
            x[pc] = 0.0;    // can be arbitrary, but set to zero
            --pc;           // , we don't have to update other things
        }
        else if(fabs(A_[pr*n+pc]) < 1e-12){
            if(fabs(x[pr]) > 1e-12){ puts("No solution"); break; }
            --pr;
        }
        else{ // update
            x[pc] /= A_[pr*n+pc];
            #pragma omp parallel for
            for(i = pr - 1; i >= 0; --i)
                x[i] -= A_[i*n+pc] * x[pc];
            --pr; --pc;
        }
    }
    // release
    free(A_); free(b_);
}


#endif
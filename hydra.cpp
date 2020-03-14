#include "hydra.h"

int all_within(int n, double *x, double lowerb, double upperb){
    int i;
    for(i = 0; i < n; ++i)
        if(isnan(x[i]) || x[i] < lowerb || x[i] > upperb)
            break;
    if(i == n)              return 1; // all are within
    else if(isnan(x[i]))    return 2; // contain NaN
    else                    return 0; // some are not within
}
void rands(int n, double *x, double a, double b){
    for(int i = 0; i < n; ++i)
        x[i] = (double)rand() / RAND_MAX * a + b;
}
void zeros(int n, double *x){
    memset(x, 0, n * sizeof(double));
}
void ones(int n, double *x){
    for(int i = 0; i < n; ++i)
        x[i] = 1.0;
}
void x_plus_sy(int n, double *x, double s, double *y, double *z){
    for(int i = 0; i < n; ++i)
        z[i] = x[i] + s * y[i];
}
void x_minus_Ay(int n, double *x, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= A[i*n+j] * y[j];
}
void x_minus_ATy(int n, double *x, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= A[j*n+i] * y[j];
}
void x_minus_sAy(int n, double *x, double s, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= s * A[i*n+j] * y[j];
}
void x_minus_sATy(int n, double *x, double s, double *A, double *y, double *z){
    memcpy(z, x, n * sizeof(double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            z[i] -= s * A[j*n+i] * y[j];
}
double xTy(int n, double *x, double *y){
    double ret = 0.0;
    for(int i = 0; i < n; ++i)
        ret += x[i] * y[i];
    return ret;
}
double xTAy(int n, double *x, double *A, double *y){
    double ret = 0.0;
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            ret += x[i] * A[i*n+j] * y[j];
    return ret;
}
void print_x(int n, double *x, const char *fmt, const char *prefix){
    if(prefix)
        printf("%s", prefix);
    for(int i = 0; i < n; ++i)
        printf(fmt, x[i]);
    putchar('\n');
}
void print_A(int n, double *A, const char *fmt, const char *prefix){
    if(prefix)
        printf("%s", prefix);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            printf(fmt, A[i*n+j]);
        putchar('\n');
    }
}

void bicg(int n, double *A, double *x, double *b, int maxattempts, int maxiters, double tol){
    int i, j, k;
    double *x_, *ri, *rj, *ri_, *rj_, *p, *p_, *tmp; // j = i + 1, _ = transpose
    double alpha, beta;

    // get parameters
    if(n * 10 > maxiters) maxiters = n * 10;

    // allocate
    x_ = (double*)malloc(n * sizeof(double));
    ri = (double*)malloc(n * sizeof(double));
    rj = (double*)malloc(n * sizeof(double));
    ri_ = (double*)malloc(n * sizeof(double));
    rj_ = (double*)malloc(n * sizeof(double));
    p = (double*)malloc(n * sizeof(double));
    p_ = (double*)malloc(n * sizeof(double));
    tmp = (double*)malloc(n * sizeof(double));

    // try different initial guess several times
    for(i = 0; i < maxattempts; ++i){
        printf("Attempt %d\n", i + 1);
        // initialize
        rands(n, x, 2.0, -1.0);
        memcpy(x_, x, n * sizeof(double));
        x_minus_Ay(n, b, A, x, ri);
        x_minus_ATy(n, b, A, x, ri_);
        memcpy(p, ri, n * sizeof(double));
        memcpy(p_, ri_, n * sizeof(double));
        // start iteration
        for(j = 0; j < maxiters; ++j){
            alpha = xTy(n, ri_, ri) / xTAy(n, p_, A, p);
            x_plus_sy(n, x, alpha, p, x);
            // check x
            x_minus_Ay(n, b, A, x, tmp); // residual
            if((k = all_within(n, tmp, -tol, tol)) > 0) break; // NaN or all within
            // update
            x_plus_sy(n, x_, alpha, p_, x_);
            x_minus_sAy(n, ri, alpha, A, p, rj);
            x_minus_sATy(n, ri_, alpha, A, p_, rj_);
            beta = xTy(n, rj_, rj) / xTy(n, ri_, ri);
            x_plus_sy(n, rj, beta, p, p);
            x_plus_sy(n, rj_, beta, p_, p_);
            memcpy(ri, rj, n * sizeof(double));
            memcpy(ri_, rj_, n * sizeof(double));
        }
        // check
        if(j == maxiters)   puts("Cannot converge.");
        else if(k == 2)     puts("Detect NaN.");
        else{               printf("Converge at iteration %d.\n", j+1); break; }
    }

    // release
    free(x_); free(ri); free(rj); free(ri_); free(rj_); free(p); free(p_); free(tmp);
}
void gaussian(int n, double *A, double *x, double *b){
    double *A_, *b_; // temp A and b to change later
    int pr, pc, pi; // pivot row, pivot col, pivot index
    int i, j, ii, jj;
    double maxabs, tmp;

    // allocate
    A_ = (double*)malloc(n * n * sizeof(double));
    b_ = (double*)malloc(n * sizeof(double));

    // initialize
    memcpy(A_, A, n * n * sizeof(double));
    memcpy(b_, b, n * sizeof(double));
    
    // To row encholon form
    for(pr = pc = 0; pr < n && pc < n; ){
        // find max abs
        pi = pr, maxabs = fabs(A_[pr*n+pc]);
        for(i = pr+1; i < n; ++i)
            if((tmp = fabs(A_[i*n+pc])) > maxabs)
                pi = i, maxabs = tmp;
        // start ellimination
        if(ISZERO(maxabs)){ // pivot is zero
            A_[pr*n+pc] = NAN; // mark as NAN, so we can recognize it later
            if(pc + 1 < n) ++pc;
            else break;
        }
        else{
            if(pi != pr){ // swap row if needed
                memcpy(x, A_+pr*n, n * sizeof(double)); // x as tmp array
                memcpy(A_+pr*n, A_+pi*n, n *sizeof(double));
                memcpy(A_+pi*n, x, n * sizeof(double));
                tmp = b_[pr], b_[pr] = b_[pi], b_[pi] = tmp;
            }
            for(i = pr + 1; i < n; ++i){ // elliminate
                tmp = -A_[i*n+pc] / A_[pr*n+pc];
                A_[i*n+pc] = 0.0;
                for(j = pc + 1; j < n; ++j)
                    A_[i*n+j] += A_[pr*n+j] * tmp;
                b_[i] += b_[pr] * tmp;
            }
            if(pc + 1 >= n || pr + 1 >= n) break;
            else{ ++pc; ++pr; }
        }
    }

    // retrieve x
    memcpy(x, b_, n * sizeof(double));
    while(pr >= 0 && pc >= 0){
        if(isnan(A_[pr*n+pc])){
            if(ISZERO(x[pc])){ // 0 x[pc] = 0
                printf("Has multiple solutions for x[%d].\n", pc);
                x[pc] = 0.0; // can be arbitrary, but set to zero to avoid updating others
            }
            else{ // 0 x[pc] = b
                puts("No solution.");
                break;
            }
            --pc;
        }
        else if(ISZERO(A_[pr*n+pc])){ // just go up if zero, then become non-zero
            --pr;
        }
        else{ // update
            x[pc] /= A_[pr*n+pc];
            for(i = pr - 1; i >= 0; --i)
                x[i] -= A_[i*n+pc] * x[pc];
            --pc;
        }
    }

    // release
    free(A_); free(b_);
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <math.h>
#include <omp.h>

// x, y, z: vector; A, B: matrix; s: scalar
inline long double swap(long double &a, long double &b){
    long double tmp = a; a = b; b = tmp;
}
inline int allzero(int n, long double *x, long double tol = 1e-12){
    int i;
    for(i = 0; i < n; ++i)
        if(isnan(x[i]) || x[i] > tol || x[i] < -tol)
            break;
    if(i == n)              return 1; // all close to 0
    else if(isnan(x[i]))    return 2; // some are nan
    else                    return 0; // some are not within
}
inline void zerox(int n, long double *x){
    memset(x, 0, n * sizeof(long double));
}
inline void randx(int n, long double *x, long double a = 1.0, long double b = 0.0){
    for(int i = 0; i < n; ++i)
        x[i] = (long double) rand() / RAND_MAX * a + b;
}
inline void randA(int n, long double *A, long double a = 1.0, long double b = 0.0){
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            A[i*n+j] = (long double) rand() / RAND_MAX * a + b;
}
inline long double sx(int n, long double s, long double *x, long double *res){
    for(int i = 0; i < n; ++i)
        res[i] = s * x[i];
}
inline long double xTy(int n, long double *x, long double *y){
    long double res = 0.;
    for(int i = 0; i < n; ++i)
        res += x[i] * y[i];
    return res;
}
inline void Ax(int n, long double *A, long double *x, long double *res){
    memset(res, 0, n * sizeof(long double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] += A[i*n+j] * x[j];
}
inline void xplusy(int n, long double *x, long double *y, long double *res){
    // #pragma omp parallel for
    for(int i = 0; i < n; ++i)
        res[i] = x[i] + y[i];
}
inline void xminusy(int n, long double *x, long double *y, long double *res){
    // #pragma omp parallel for
    for(int i = 0; i < n; ++i)
        res[i] = x[i] - y[i];
}
inline void xplussy(int n, long double *x, long double s, long double *y, long double *res){
    for(int i = 0; i < n; ++i)
        res[i] = x[i] + s * y[i];
}
inline void axplusby(int n, long double a, long double *x, long double b, long double *y, long double *res){
    // #pragma omp parallel for
    for(int i = 0; i < n; ++i)
        res[i] = a * x[i] + b * y[i];
}
inline void xminussy(int n, long double *x, long double s, long double *y, long double *res){
    // #pragma omp parallel for
    for(int i = 0; i < n; ++i)
        res[i] = x[i] - s * y[i];
}
inline void xminusAy(int n, long double *x, long double *A, long double *y, long double *res){
    memcpy(res, x, n * sizeof(long double));
    for(int i = 0; i < n; ++i)
        for(int j = 0; j < n; ++j)
            res[i] -= A[i*n+j] * y[j];
}
void printx(const char *prefix, int n, long double *x){
    printf("%s", prefix);
    for(int i = 0; i < n; ++i)
        printf("%.8Lf ", x[i]);
    putchar('\n');
}
void printA(const char *prefix, int n, long double *A){
    printf("%s", prefix);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            printf("%.3Lf ", A[i*n+j]);
        putchar('\n');
    }
}
void loadx(const char *filename, int &n, long double *&x){
    FILE *fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "Failed to load %s\n", filename);
        return;
    }
    fscanf(fp, "%d", &n);
    x = (long double*) malloc(n * sizeof(long double));
    for(int i = 0; i < n; ++i)
        fscanf(fp, "%Lf", &x[i]);
    fclose(fp);
}
void loadA(const char *filename, int &n, long double *&A){
    FILE *fp = fopen(filename, "w");
    if(fp == NULL){
        fprintf(stderr, "Failed to load %s\n", filename);
        return;
    }
    fscanf(fp, "%d", &n);
    A = (long double*) malloc(n * n * sizeof(long double));
    for(int n2 = n * n, i = 0; i < n2; ++i)
        fscanf(fp, "%Lf", &A[i]);
    fclose(fp);
}
void exportx(const char *filename, int n, long double *x){
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d\n\n", n);
    for(int i = 0; i < n; ++i)
        fprintf(fp, "%Lf\n", x[i]);
    fclose(fp);
}
void exportA(const char *filename, int n, long double *A){
    FILE *fp = fopen(filename, "w");
    fprintf(fp, "%d\n\n", n);
    for(int i = 0; i < n; ++i){
        for(int j = 0; j < n; ++j)
            fprintf(fp, "%Lf ", A[i*n+j]);
        fputc('\n', fp);
    }
    fclose(fp);
}

// Noets: for large matrix, avoid using small tol
//        Otherwise, it would be hard to converge
void bicgstab(int n, long double *A, long double *b, long double *x, bool precond = false, long double tol = 1e-6, int maxiter = 10000, int maxattempt = 1000){
    int i, j, k, pr, pc, pi, ii, jj;
    long double *A0, *b0, *r0, *r, *p, *v, *s, *t, *tmp;
    long double gmax, dtmp, rhoi, rhoj, alpha, beta, omega; // j = i + 1
    // adjust parameters
    if(10 * n > maxiter)     maxiter = 10 * n;
    if(100 * n > maxattempt) maxattempt = 100 * n;
    // allocate memory
    A0 = precond? (long double*)malloc(n * n * sizeof(long double)) : A;
    b0 = precond? (long double*)malloc(n * sizeof(long double)) : b;
    r0 = (long double*) malloc(n * sizeof(long double));
    r = (long double*) malloc(n * sizeof(long double));
    p = (long double*) malloc(n * sizeof(long double));
    v = (long double*) malloc(n * sizeof(long double));
    s = (long double*) malloc(n * sizeof(long double));
    t = (long double*) malloc(n * sizeof(long double));
    tmp = (long double*) malloc(n * sizeof(long double));
    // preconditioning
    if(precond){
        memcpy(A0, A, n * n * sizeof(long double));
        memcpy(b0, b, n * sizeof(long double));
        for(pr = pc = 0; pr < n && pc < n; ){
            // find max abs
            pi = pr;
            gmax = fabs(A0[pr*n+pc]);
            for(i = pr+1; i < n; ++i){
                dtmp = fabs(A0[i*n+pc]);
                if(dtmp > gmax) pi = i, gmax = dtmp;
            }
            // start ellimination
            if(gmax < 1e-16){ // pivot is zero
                A0[pr*n+pc] = NAN; // for convenience later
                ++pc;
            }
            else{
                // swap row
                if(pr != pi){
                    for(ii = pr * n, jj = pi * n, i = pc; i < n; ++i)
                        swap(A0[ii+i], A0[jj+i]);
                    swap(b0[pr], b0[pi]);
                }
                // elliminate
                for(i = pr + 1; i < n; ++i){
                    dtmp = A0[i*n+pc] / A0[pr*n+pc];
                    A0[i*n+pc] = 0.0;
                    for(j = pc + 1; j < n; ++j)
                        A0[i*n+j] -= A0[pr*n+j] * dtmp;
                    b0[i] -= b0[pr] * dtmp;
                }
                ++pc;
                ++pr;               
            }
        }
    }
    // several attempts with different initial guess
    for(i = 0; i < maxattempt; ++i){
        printf("Attempt %d\n", i + 1);
        // init
        if(i == 0)                       zerox(n, x);
        else if(i * 3 < maxattempt)      randx(n, x, 2.0, -1.0);
        else if(i * 3 < maxattempt * 2)  randx(n, x, 200.0, -100.0);
        else                             randx(n, x, 0.02, -0.01);
        xminusAy(n, b0, A0, x, r);
        memcpy(r0, r, n * sizeof(long double));
        // dtmp = 0.1 * sqrt(xTy(n, r, r));
        // randx(n, tmp, dtmp, -0.5 * dtmp);
        // xplusy(n, r, tmp, r0);
        rhoi = alpha = omega = 1.0;
        zerox(n, p);
        zerox(n, v);
        // start iteration
        for(j = 0; j < maxiter; ++j){
            // stage 1
            rhoj = xTy(n, r0, r);
            if((rhoj < 1e-12 && rhoj > -1e-12) || (rhoj < -1e12 || rhoj > 1e12))
                rhoj = rhoi; // avoid rho become too small or too large
            beta = rhoj * alpha / rhoi / omega;
            xminussy(n, p, omega, v, tmp);
            xplussy(n, r, beta, tmp, p);
            Ax(n, A0, p, v);
            alpha = rhoj / xTy(n, r0, v);
            xplussy(n, x, alpha, p, x);
            // check
            xminusAy(n, b0, A0, x, tmp);
            if((k = allzero(n, tmp, tol)) > 0) break;
            // stage2
            xminussy(n, r, alpha, v, s);
            Ax(n, A0, s, t);
            omega = xTy(n, t, x) / xTy(n, t, t);
            xplussy(n, x, omega, s, x);
            // check
            xminusAy(n, b0, A0, x, tmp);
            if((k = allzero(n, tmp, tol)) > 0) break;
            // update
            xminussy(n, s, omega, t, r);
            rhoi = rhoj;
        }
        // check
        if(j == maxiter)
            puts("Caonnot converge");
        else if(k == 2) // break due to nan
            puts("Detect NaN");
        else{ // converge
            printf("Converge at iteration %d\n", j + 1);
            break;
        }
    }

    // release
    if(precond){ free(A0); free(b0); }
    free(r0); free(r); free(p); free(v); free(s); free(t); free(tmp);
}

void gaussian(int n, long double *A, long double *b, long double *x){
    long double *G;
    int pr, pc, pi; // pivot row, pivot column, pivot index
    int i, j, ii, jj;
    long double gmax, tmp;

    // allocate
    G = (long double*) malloc(n * n * sizeof(long double));
    // init
    memcpy(G, A, n * n * sizeof(long double));
    memcpy(x, b, n * sizeof(long double)); // x as temporary b
    // To row encholon form
    for(pr = pc = 0; pr < n && pc < n; ++pc){
        // find max abs
        pi = pr, gmax = fabs(G[pr*n+pc]);
        for(i = pr+1; i < n; ++i)
            if((tmp = fabs(G[i*n+pc])) > gmax)
                pi = i, gmax = tmp;
        // start ellimination
        if(gmax < 1e-12)        // pivot is zeri
            G[pr*n+pc] = NAN;   // for convenience later
        else{
            if(pr != pi){ // swap row
                ii = pr*n, jj = pi*n;
                // #pragma omp parallel for
                for(i = pc; i < n; ++i)
                    swap(G[ii+i], G[jj+i]);
                swap(x[pr], x[pi]);
            }
            for(i = pr + 1; i < n; ++i){ // elliminate
                tmp = G[i*n+pc] / G[pr*n+pc];
                G[i*n+pc] = 0.0;
                // #pragma omp parallel for
                for(j = pc + 1; j < n; ++j)
                    G[i*n+j] -= G[pr*n+j] * tmp;
                x[i] -= x[pr] * tmp;
            }
            ++pr;
        }
    }
    // retrieve x
    for(pr = pc = n - 1; pr >= 0 && pc >=0; ){
        if(pr + 1 < n && isnan(G[(pr+1)*n+pc])){
            printf("Has multiple solutions for x[%d]\n", pc);
            x[pc] = 0.0;    // can be arbitrary, but set to zero
            --pc;           // , we don't have to update other things
        }
        else if(fabs(G[pr*n+pc]) < 1e-12){
            if(fabs(x[pr]) > 1e-12){ puts("No solution"); break; }
            --pr;
        }
        else{ // update
            x[pc] /= G[pr*n+pc];
            // #pragma omp parallel for
            for(i = pr - 1; i >= 0; --i)
                x[i] -= G[i*n+pc] * x[pc];
            --pr; --pc;
        }
    }
    // release
    free(G);
}

void gaussian2(int n, long double *A, long double *b, long double *x, long double *bufferx){
    long double *G;
    int pr, pc, pi; // pivot row, pivot column, pivot index
    int i, j, ii, jj;
    long double gmax, tmp;

    // allocate
    G = (long double*) malloc(n * n * sizeof(long double));
    // init
    memcpy(G, A, n * n * sizeof(long double));
    memcpy(x, b, n * sizeof(long double)); // x as temporary b
    // To row encholon form
    for(pr = pc = 0; pr < n && pc < n; ++pc){
        // find max abs
        pi = pr, gmax = fabs(G[pr*n+pc]);
        for(i = pr+1; i < n; ++i)
            if((tmp = fabs(G[i*n+pc])) > gmax)
                pi = i, gmax = tmp;
        // start ellimination
        if(gmax < 1e-12)        // pivot is zeri
            G[pr*n+pc] = NAN;  // for convenience later
        else{
            if(pi != pr){ // swap row
                memcpy(bufferx, G+pr*n, n * sizeof(long double));
                memcpy(G+pr*n, G+pi*n, n * sizeof(long double));
                memcpy(G+pi*n, bufferx, n * sizeof(long double));
                swap(x[pr], x[pi]);
            }
            for(i = pr + 1; i < n; ++i){ // elliminate
                tmp = -G[i*n+pc] / G[pr*n+pc];
                G[i*n+pc] = 0.0;
                #pragma omp parallel for
                for(j = pc + 1; j < n; ++j)
                    G[i*n+j] += G[pr*n+j] * tmp;
                x[i] += x[pr] * tmp;
            }
            ++pr;
        }
    }
    // retrieve x
    for(pr = pc = n - 1; pr >= 0 && pc >=0; ){
        if(pr + 1 < n && isnan(G[(pr+1)*n+pc])){
            printf("Has multiple solutions for x[%d]\n", pc);
            x[pc] = 0.0;    // can be arbitrary, but set to zero
            --pc;           // , we don't have to update other things
        }
        else if(fabs(G[pr*n+pc]) < 1e-12){
            if(fabs(x[pr]) > 1e-12){ puts("No solution"); break; }
            --pr;
        }
        else{ // update
            x[pc] /= G[pr*n+pc];
            #pragma omp parallel for
            for(i = pr - 1; i >= 0; --i)
                x[i] -= G[i*n+pc] * x[pc];
            --pr; --pc;
        }
    }
    // release
    free(G);
}
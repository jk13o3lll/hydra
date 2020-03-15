#include "hydra_cuda.h"
#include <curand.h>

// https://docs.nvidia.com/cuda/cuda-math-api/index.html
// https://docs.nvidia.com/cuda/cusolver/index.html

// gpu kernels for bicg
void all_within(int n, double *x, double lowerb, double upperb){
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    // do this outside? no good way to do this in GPU?
}
__global__ void rands(int n, double *x, double a = 1.0, double b = 0.0){
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    // https://stackoverflow.com/questions/22425283/how-could-we-generate-random-numbers-in-cuda-c-with-different-seed-on-each-run
    // https://developer.nvidia.com/curand
}
__global__ void x_plus_sy(int n, double *x, double s, double *y, double *z){ // minus by -s
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id < n)
        z[id] = x[id] + s * y[id];
}
__global__ void Ax(int n, double *A, double *x, double *y){ // y = Ax
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id < n){
        double tmp = 0.0;
        for(int i = 0; i < n; ++i)
            tmp += A[id*n+i] * x[i];
        y[id] = tmp;
    }
}
__global__ void ATx(int n, double *A, double *x, double *y){ // y = ATx
    int id = blockDim.x * blockIdx.x + threadIdx.x;
    if(id < n){
        double tmp = 0.0;
        for(int i = 0; i < n; ++i)
            tmp += A[i*n+id] * x[i];
        y[id] = tmp;
    }
}
__global__ void xTy(int n, float *x, float *y, float *z){ // z = xTy
    __shared__ float tmp[TPB];  // shared memory useful when data is used serveral times
    int id = blockDim.x*blockIdx.x + threadIdx.x;
    int ii = threadIdx.x;

    tmp[ii] = id<n? x[id]*y[id] : 0.0;
    z[id] = 0.0;
    __syncthreads();
    for(int stride = blockDim.x >> 1; stride > 0; stride >>= 1){ // parallel reduction
        if(ii < stride) tmp[ii] += tmp[ii + stride]; // blockDim.x should be 2^n, otherwise some entries will miss
        __syncthreads();
    }
    if(ii == 0) z[blockIdx.x] = tmp[0];
    __syncthreads();
    if(id == 0){
        float sum = 0.0;
        for(int i = 0; i < gridDim.x; i += 2) sum += z[i];
        z[0] = sum;
    }
    else if(id == 1){
        float sum = 0.0;
        for(int i = 1; i < gridDim.x; i += 2) sum += z[i];
        z[1] = sum;
    }
    __syncthreads();
    if(id == 0) z[0] += z[1];
}
void bicg(int n, double *A, double *x, double *b, int maxattempts = 100, int maxiters = 1000, double tol = 1e-6){
    int i, j, k;

    if(n * 10 > maxiters) maxiters = n * 10;
    // malloc in solve_cuda
    // try different initial guess several times
    for(i = 0; i < maxattemps; ++i){
        // init
        // rands(n, x, n * sizeof(double));
        // ...
        // start iteration
        for(j = 0; j < maxiters; ++j){
            // ...

            // don't check residual every iteration
            // for example, check every 100 iteration

            // ...
        }
        // check

    }
    // free in solve_cuda
}


__global__ void compute_R_and_J(){
    int id = blockDim.x*blockIdx.x + threadIdx.x;
    // compute R

    // compute J
}
void solve_cuda(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiters, int maxattempts, double tol, double step){
    int i, j;

    // allocate memory

    // send to device

    // newton's method
    for(i = 0; i < maxattempts; ++i){
        // init cpu vaiables

        // init gpu variables

        // start iterations
        for(j = 0; j < maxiters; ++j){


        }
        // check
    }

    // get from device

    // free_memory

}

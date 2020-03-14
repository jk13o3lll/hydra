#ifndef __HYDRA_CUDA_H__
#define __HYDRA_CUDA_H__

// #define N 200
#define TPB 32 // threads per block
// #define BPG ((N+TPB-1)/TPB) // blocks per grid

void solve_cuda(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiters, int maxattempts, double tol, double step);

#endif
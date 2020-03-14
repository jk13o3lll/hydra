#ifndef __HYDRA_CUDA_H__
#define __HYDRA_CUDA_H__









// compute_R_and_J and bicg are parallelized by CUDA
void solve_cuda(int nE, int nLeq, int nNeq, double *incLoop, double *conLoop, double *incNode, double *conNode, double *&x, double n, int maxiter, int maxattempts, double tol, double step);

#endif
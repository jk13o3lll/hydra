#include "cuda_main.h"

// === Pipe flow ===
float *h_rad, *h_Q, *h_Rcc;
float *d_rad, *d_Q, *d_Rcc;
int *h_Rc, *h_nRc;
int *d_Rc, *d_nRc;
float *d_R, *d_J, *d_negdQ, *d_residual_NR;

// === CG ===
float *d_x, *d_r, *d_p;
float *d_x_, *d_r_, *d_p_;
float *d_scalars; // p_Ap, r_r old, r_r new, alpha, beta
float *d_Ap, *d_p_A;

void allocate_memory(){
    size_t sz = N * sizeof(float);
    cudaError_t err;
    h_a = (float*) malloc(sz);
    err = cudaMalloc((void**)&d_a, sz);
    printf("CUDA error (malloc d_a) = %s\n", cudaGetErrorString(err));
}
void free_memory(){
    if(h_a) free(h_a);
    if(d_a) cudaFree(d_a);
}

void input_from_file(){
    for(int i=0; i<N; ++i)
        h_a[i] = 1.0;
        // h_a[i] = i;
}
void output_to_file(){
    for(int i=0; i<N; ++i)
        printf("%f ", h_a[i]);
    putchar('\n');
    // printf("%f\n", h_a[0]);
}

void send_to_device(){
    size_t sz = N * sizeof(float);
    cudaError_t err;

    err = cudaMemcpy(d_a, h_a, sz, cudaMemcpyHostToDevice);
    printf("CUDA error (memcpy h_a -> d_a) = %s\n", cudaGetErrorString(err));
}
void get_from_device(){
    size_t sz = N * sizeof(float);
    cudaError_t err;

    err = cudaMemcpy(h_a, d_a, sz, cudaMemcpyDeviceToHost);
    printf("CUDA error (memcpy d_a -> h_a) = %s\n", cudaGetErrorString(err));
}


__global__ void DotProduct(float *x, float *y, float *z, int n){ // z = xTy
    __shared__ float tmp[TPB];  // shared memory useful when data is used serveral times
    int id = blockDim.x*blockIdx.x + threadIdx.x;
    int ii=threadIdx.x;

    tmp[ii] = id<n? x[id]*y[id] : 0.0;
    z[id] = 0.0;
    __syncthreads();
    for(int stride = blockDim.x >> 1; stride > 0; stride >>= 1){ // parallel reduction
        if(ii<stride) tmp[ii] += tmp[ii + stride];             // blockDim.x should be 2^n, otherwise some entry will miss
        __syncthreads();
    }
    if(ii == 0) atomicAdd(&z[0], tmp[0]);
}
__global__ void MatrixVectorProduct(float *A, float *x, float *b, int n){ // Ax=b
	int id = blockDim.x*blockIdx.x + threadIdx.x;

	if(id < n){
		float sum = 0.0;
		for(int i=0; i<n; ++i)	sum += A[id*n+i] * x[i];
		b[id] = sum;
	}
}
__global__ void VectorMatrixProduct(float *x, float *A, float *b, int n){ // ATx=b
	int id = blockDim.x*blockIdx.x + threadIdx.x;

	if(id<n){
		float sum = 0.0;
		for(int i=0; i<n; ++i)	sum += A[i*n+id] * x[i];
		b[id] = sum;
	}
}

__global__ void BiCG_Init(float *b, float *x, float *x_, float *r, float *r_, float *p, float *p_, int n){
	// set both x to zero
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	if(id<n){
		x[id] = 0;	x_[id] = 0;
		r[id] = b[id];	r_[id] = b[id];
		p[id] = b[id];	p_[id] = b[id];
	}
}
__global__ void BiCG_Compute_Alpha(float *scalars){
	scalars[3] = scalars[1]/scalars[0];
}
__global__ void BiCG_Update_x_and_r(float *x, float *x_, float *p, float *p_, float *r, float *r_, float *Ap, float *p_A, float *scalars, int n){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	if(id<n){
		float alpha = scalars[3];
		x[id] += alpha*p[id];	x_[id] += alpha*p_[id];
		r[id] -= alpha*Ap[id];	r_[id] -= alpha*p_A[id];
	}
}
__global__ void BiCG_Compute_Beta(float *scalars){
	scalars[4] = scalars[2]/scalars[1];
}
__global__ void BiCG_Update_p(float *p, float *p_, float *r, float *r_, float *scalars, int n){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	if(id<n){
		float beta = scalars[4];
		p[id] = r[id] + beta*p[id];
		p_[id] = r_[id] + beta*p_[id];
	}
}
__global__ void BiCG_Update_residual(float *scalars){
	scalars[1] = scalars[2];
}
void BiCG(float *A, float *x, float *b, int n){
	int i;
	float residual; // RTR

	BiCG_Init<<<BPG,TPB>>>(b, x, d_x_, d_r, d_r_, d_p, d_p_, n);
	DotProduct<<<BPG,TPB>>>(d_r, d_r, d_scalars+1, n, d_tmp_BPG);
	for(i=0; i<ITER_CG; ++i){
		MatrixVectorProduct<<<BPG,TPB>>>(A, d_p, d_Ap, n);
		VectorMatrixProduct<<<BPG,TPB>>>(d_p_, A, d_p_A, n);
		DotProduct<<<BPG,TPB>>>(d_p_, d_Ap, d_scalars, n, d_tmp_BPG);
		BiCG_Compute_Alpha<<<1,1>>>(d_scalars);
		BiCG_Update_x_and_r<<<BPG,TPB>>>(x, d_x_, d_p, d_p_, d_r, d_r_, d_Ap, d_p_A, d_scalars, n);
		DotProduct<<<BPG,TPB>>>(d_r_, d_r, d_scalars+2, n, d_tmp_BPG); // RTR new
		cudaMemcpy(&residual, d_scalars+2, sizeof(float), cudaMemcpyDeviceToHost);
		if(fabsf(residual) < EPS_CG)	break;
		BiCG_Compute_Beta<<<1,1>>>(d_scalars);
		BiCG_Update_p<<<BPG,TPB>>>(d_p, d_p_, d_r, d_r_, d_scalars, n);
		BiCG_Update_residual<<<1,1>>>(d_scalars);
	}
	//printf("i=%d\n", i);
}

// Pipe flow
__global__ void Compute_R_and_J(float *rad, float *Q, float *Rcc, int *Rc, int *nRc, float *R, float *J, int n){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	int i, tmp, offset = id*n;
	float sum, tmpQ;

	if(id < N_LOOPS){ // loop eqns
		sum = Rcc[id];
		for(i=0; i<n; ++i)	J[offset+i] = 0;
		for(i=0; i<nRc[id]; ++i){
			tmp = Rc[offset+i];
			if(tmp < 0){
				tmp = -tmp-1;
				tmpQ = Q[tmp];
				sum -= rad[tmp]*fabsf(tmpQ)*tmpQ; // R
				J[offset+tmp] = -2*rad[tmp]*tmpQ; // J
			}
			else{
				tmp = tmp-1;
				tmpQ = Q[tmp];
				sum += rad[tmp]*fabsf(tmpQ)*tmpQ;
				J[offset+tmp] = 2*rad[tmp]*tmpQ;
			}
		}
		R[id] = sum;
	}
	else if(id < n){ // nodes eqns
		sum = Rcc[id];
		for(i=0; i<n; ++i)	J[offset+i] = 0;
		for(i=0; i<nRc[id]; ++i){
			tmp = Rc[offset+i];
			if(tmp < 0){
				tmp = -tmp-1;
				sum -= Q[tmp]; // R
				J[offset+tmp] = -1.0; // J
			}
			else{
				tmp = tmp-1;
				sum += Q[tmp];
				J[offset+tmp] = 1.0;
			}
		}
		R[id] = sum;
	}
}
__global__ void Update_Q(float *Q, float *negdQ, int n){
	int id = blockDim.x*blockIdx.x + threadIdx.x;
	if(id<n){
		Q[id] -= DAMP_NR*negdQ[id];
	}
}

void Newton_Raphson(float *R, float *J, float *Q, int n){ // find the sol of R=0, with NR use J*dx=-R, x=x+dx => J*(-dx)=R, x=x-(-dx)
	int i;
	float residual;

	for(i=0; i<ITER_NR; ++i){
		Compute_R_and_J<<<BPG,TPB>>>(d_rad, Q, d_Rcc, d_Rc, d_nRc, R, J, n);
		DotProduct<<<BPG,TPB>>>(R, R, d_residual_NR, n, d_tmp_BPG);
		cudaMemcpy(&residual, d_tmp_BPG, sizeof(float), cudaMemcpyDeviceToHost);
		if(residual < EPS_NR)	break;
		BiCG(J, d_negdQ, R, n);
		Update_Q<<<BPG,TPB>>>(Q, d_negdQ, n);
	}

	printf("iteration_NR = %d\n", i);
	printf("residual_NR = %.10e\n", residual);
}

// https://www.cs.cmu.edu/afs/cs/academic/class/15668-s11/www/cuda-doc/html/group__CUDART__DEVICE_g5aa4f47938af8276f08074d09b7d520c.html
// https://zhuanlan.zhihu.com/p/41151532



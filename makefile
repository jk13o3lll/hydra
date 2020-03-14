CUDA_INC_DIR = /usr/local/cuda/include # nvcc knows the path
CUDA_LIB_DIR = /usr/local/cuda/lib64
CUDA_LIBS = -lcuda -lcudart -lcurand # -lcusparse -lcusolver

# make CPU version
cpu:
	g++ -O3 main.cpp hydra.cpp -o test_cpu.exe

# make CUDA version
cuda: cpu_core gpu_core
	g++ main.o hydra.o hydra_cuda.o -L$(CUDA_LIB_DIR) $(CUDA_LIBS) -O3 -o test_cuda.exe

cpu_core:
	g++ -O3 -c hydra.cpp -D HYDRA_USE_CUDA
	g++ -O3 -c main.cpp -D HYDRA_USE_CUDA

gpu_core:
	nvcc -O3 -c hydra_cuda.cu

# clean
.PHONY: clean
clean:
	rm *.o

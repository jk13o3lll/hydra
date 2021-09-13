# CUDA_INC_DIR = /usr/local/cuda/include # cvaa knows the path
CUDA_LIB_DIR = /usr/local/cuda/lib64
CUDA_LIBS = -lcuda -lcudart -lcusparse -lcusolver

all: MAIN

MAIN: CPU_CORE GPU_CORE
	g++ main.o cuda_main.o -L$(CUDA_LIB_DIR) $(CUDA_LIBS) -O3 -o test.exe

CPU_CORE:
	g++ -O3 -c main.c

GPU_CORE:
	nvcc -O3 -c cuda_main.cu

clean:
	rm *.o

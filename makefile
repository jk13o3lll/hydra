CUDA_INC_DIR = /usr/local/cuda/include
CUDA_LIB_DIR = /usr/local/cuda/lib64
CUDA_LIBS = -lcudart -lcusparse -lcusolver

CPU: # cpu version
# g++ compile all

GPU: # cuda version
# nvcc compile cuda files first
# g++ compile main function and link all

clean:

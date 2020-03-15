#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define N 200
#define TPB 32 // threads per block
#define BPG ((N+TPB-1)/TPB) // blocks per grid

void allocate_memory();
void free_memory();

void input_from_file();
void output_to_file();

void send_to_device();
void get_from_device();

void call_gpu_kernels();

#include "cuda_main.h"

int main(){
    allocate_memory();
    input_from_file();
    send_to_device();
    call_gpu_kernels();
    get_from_device();
    output_to_file();
    free_memory();
    return 0;
}

// Allocate_Memory();
// Input_From_File();
// Send_To_Device();
// Call_GPU_Kernals();
// Get_From_Device();
// Output_To_File();
// Free_Memory();
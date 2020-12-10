#include "CudaRNG.h"
#include <cuda.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include <curand.h>
#include <curand_kernel.h>
#include <Windows.h>
#include <omp.h>

//http://ianfinlayson.net/class/cpsc425/notes/cuda-random

__global__ void init(unsigned int seed, curandState_t* states) {
	curand_init(seed,
		blockIdx.x,
		0,
		&states[blockIdx.x]);
}

/* this GPU kernel takes an array of states, and an array of ints, and puts a random int into each */
__global__ void randoms(curandState_t* states, unsigned int* numbers) {
	numbers[blockIdx.x] = curand(&states[blockIdx.x]);
}


void CudaRNG::getRandomNums(const int bufferSize, int * buffer) {
	#pragma omp critical GPUAccess
	{
		curandState_t* states;
		cudaMalloc((void**)&states, bufferSize * sizeof(curandState_t));
		init <<<bufferSize, 1 >> > (GetTickCount(), states);

		unsigned int* gpu_nums;
		cudaMalloc((void**)&gpu_nums, bufferSize * sizeof(unsigned int));

		randoms <<<bufferSize, 1 >> > (states, gpu_nums);
		cudaMemcpy(buffer, gpu_nums, bufferSize * sizeof(unsigned int), cudaMemcpyDeviceToHost);

		cudaFree(states);
		cudaFree(gpu_nums);
	}
}
#include "BufferedRandomNumberGenerator.h"

#include <cuda.h>
#include <curand.h>
#include <curand_kernel.h>
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

__global__ void generateRandomInts(unsigned int seed, int* result) {
	curandState_t state;
	curand_init(seed, /* the seed controls the sequence of random values that are produced */
		blockIdx.x, /* the sequence number is only important with multiple cores */
		0, /* the offset is how much extra we advance in the sequence for each call, can be 0 */
		&state);
	*result = curand(&state);

}
/*
void CudaRNG::fillBuffer(int * buffer, int bufferSize) {
	//BufferedRandomNumberGenerator::mtx.lock();

	int* gpu_x;
	cudaMalloc((void**)&gpu_x, sizeof(int));
	//generateRandomInts<<<1, 1 >>>(123, gpu_x);

	int x;
	cudaMemcpy(&x, gpu_x, sizeof(int), cudaMemcpyDeviceToHost);
	printf("Random number = %d.\n", x);
	cudaFree(gpu_x);


	//BufferedRandomNumberGenerator::mtx.unlock();
}*/
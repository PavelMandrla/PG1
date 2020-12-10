#pragma once

#include "CudaRNG.h"


class RNG {
private:
	int currI;

	int bufferSize;
	int *buffer;
	std::shared_ptr<CudaRNG> cudaRNG;

	std::mt19937 generator;
public:
	RNG(int bufferSize, std::shared_ptr<CudaRNG> cudaRNG);
	~RNG();
	float getRandNum(float from, float to);

};


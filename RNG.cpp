#include "stdafx.h"
#include "RNG.h"
#include <Windows.h>

RNG::RNG(int bufferSize, std::shared_ptr<CudaRNG> cudaRNG) {
	this->bufferSize = bufferSize;
	this->cudaRNG = cudaRNG;
	this->buffer = new int[this->bufferSize];
	
	//this->cudaRNG->getRandomNums(this->bufferSize, this->buffer);
	for (int i = 0; i < bufferSize; i++) {
		this->buffer[i] = rand();
	}
	this->currI = 0;
	this->generator = std::mt19937(GetTickCount());
}

RNG::~RNG() {
	delete this->buffer;
}

float RNG::getRandNum(float from, float to) {
	std::uniform_real_distribution<float> uni_dist(from, to);
	return uni_dist(this->generator);


	/*float result = 1.0f;
	if (this->currI < this->bufferSize) {
		//result = float(this->buffer[currI]) / 1000.0f;
		result = from + float(this->buffer[currI]) / ((float(INT_MAX) / (to - from)));
	}
	this->currI++;
	return result;*/
	//return 1.0f;
	/*
	float result = from + static_cast <float> (this->buffer[currI]) / (static_cast <float> (INT_MAX / (to - from)));
	this->currI++;
	if (this->currI >= this->bufferSize) {
		//this->cudaRNG->getRandomNums(this->bufferSize, this->buffer);
		this->currI = 0;
	}
	return result;
	*/
}

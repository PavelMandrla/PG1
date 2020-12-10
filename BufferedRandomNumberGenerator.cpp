#include "stdafx.h"
#include "BufferedRandomNumberGenerator.h"



BufferedRandomNumberGenerator::BufferedRandomNumberGenerator(int bufferSize) {
	this->bufferSize = bufferSize;
	this->buffer = new int[this->bufferSize];
	//BufferedRandomNumberGenerator::fillBuffer(this->buffer, this->bufferSize);
	this->currentIndex = 0;
}

BufferedRandomNumberGenerator::~BufferedRandomNumberGenerator() {
	delete this->buffer;
}

float BufferedRandomNumberGenerator::getRandFloat(float from, float to) {
	//BufferedRandomNumberGenerator::mtx.lock();
	float range = (to - from);
	float div = RAND_MAX / range;
	float randNum =  from + (this->buffer[currentIndex] / div);

	currentIndex++;
	if (currentIndex = this->bufferSize) {
		//BufferedRandomNumberGenerator::fillBuffer(this->buffer, this->bufferSize);
		this->currentIndex = 0;
	}
	
	return randNum;
}



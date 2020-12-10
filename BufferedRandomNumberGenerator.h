#pragma once
#include <memory>

class BufferedRandomNumberGenerator {
private:
	int currentIndex;
	int bufferSize;
	int *buffer;
public:
	//static std::mutex mtx;
	BufferedRandomNumberGenerator(int bufferSize);
	~BufferedRandomNumberGenerator();
	float getRandFloat(float from, float to);

};

class Singleton {
private:
	/* Here will be the instance stored. */
	static Singleton* instance;

	/* Private constructor to prevent instancing. */
	Singleton();

public:
	/* Static access method. */
	static Singleton* getInstance();
};

Singleton* Singleton::instance = 0;

Singleton* Singleton::getInstance() {
	if (instance == 0) {
		instance = new Singleton();
	}

	return instance;
}

Singleton::Singleton() {
}
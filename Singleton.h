#pragma once
class Singleton {
	private:
		/* Here will be the instance stored. */
		static Singleton* instance;

		/* Private constructor to prevent instancing. */
		Singleton();

	public:
		/* Static access method. */
		static Singleton* getInstance();

	/* Null, because instance will be initialized on demand. */
	
};


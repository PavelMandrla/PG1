#include "stdafx.h"
#include "tutorials.h"

int main()
{
	printf( "PG1, (c)2011-2019 Tomas Fabian\n\n" );

	_MM_SET_FLUSH_ZERO_MODE( _MM_FLUSH_ZERO_ON );
	_MM_SET_DENORMALS_ZERO_MODE( _MM_DENORMALS_ZERO_ON );

	//return tutorial_1();
	//return tutorial_2();

	return raytracerDemo( "../../../data/6887_allied_avenger.obj" );
	//return pathtracerDemo("../../../data/cornell_box2.obj");
}

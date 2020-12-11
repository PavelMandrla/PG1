#pragma once

//#include "simpleguidx11.h"
#include "Tracer.h"
#include "corecrt_math_defines.h"
#include "vector3.h"
#include "structs.h"
#include "surface.h"
#include "camera.h"

#include <string>


class PathTracer : public Tracer {
public:
	PathTracer(const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=4,verbose=3");

	Color4f get_pixel(const int x, const int y, const float t) override;
private:
	void getCosWeightedSample(Vector3 n, Vector3 &omega_i, float &pdf);
	Vector3 rotateVector(Vector3 v, Vector3 n);

	

	virtual Color4f trace(RTCRay ray, int level, float rayIOR = 1) override;
};


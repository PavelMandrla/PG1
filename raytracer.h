#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "LightSource.h"
#include "SphericalMap.h"
#include "BVH.h"
#include "corecrt_math_defines.h"
#include "Tracer.h"
#include <omp.h>
#include <random>
#include <array>

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public Tracer {
public:
	Raytracer( const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config = "threads=0,verbose=3" );

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

private:
	LightSource* ambient;
	std::vector<LightSource> lightSources;

	std::shared_ptr<BVH> bvh;
	std::shared_ptr<SphericalMap> sphericalMap;

	bool isBlocked(LightSource lightSrc, Vector3 hitPoint);
	Color4f getPhongIllumination(Vector3 hitPoint, Vector3 d, Vector3 n, Material *material, RTCGeometry geometry, RTCRayHit ray_hit);

	virtual Color4f trace(RTCRay ray, int level, float rayIOR = 1) override;
};

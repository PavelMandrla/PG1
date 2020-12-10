#pragma once
#include "simpleguidx11.h"
#include "surface.h"
#include "camera.h"
#include "LightSource.h"
#include "SphericalMap.h"
#include "BVH.h"
#include "corecrt_math_defines.h"
#include <omp.h>
#include <random>
#include <array>

/*! \class Raytracer
\brief General ray tracer class.

\author Tomáš Fabián
\version 0.1
\date 2018
*/
class Raytracer : public SimpleGuiDX11
{
public:
	Raytracer( const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3" );
	~Raytracer();

	int InitDeviceAndScene( const char * config );

	int ReleaseDeviceAndScene();

	void LoadScene( const std::string file_name );

	Color4f get_pixel( const int x, const int y, const float t = 0.0f ) override;

	int Ui();

private:
	LightSource* ambient;
	std::vector<LightSource> lightSources;
	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	std::shared_ptr<BVH> bvh;

	std::shared_ptr<SphericalMap> sphericalMap;

	Color4f trace(RTCRay ray, int level, float rayIOR = 1);
	Color4f pathTrace(RTCRay ray, int level, float rayIOR = 1);
	bool isBlocked(LightSource lightSrc, Vector3 hitPoint);


	RTCDevice device_;
	RTCScene scene_;
	Camera camera_;

	RTCRay getRay(Vector3 origin, Vector3 direction);
	Vector3 getRefractedVector(Vector3 d, Vector3 n, float ior1, float ior2);
	Vector3 getReflectedVector(Vector3 d, Vector3 n);

	float getRaySegmentLength(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR);
	float getReflectionCoefficient(Vector3 d, Vector3 n, float ior1, float ior2);

	Color4f getReflectedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level);
	Color4f getRefractedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level);
	
	Color4f getPhongIllumination(Vector3 hitPoint, Vector3 d, Vector3 n, Material *material, RTCGeometry geometry, RTCRayHit ray_hit);

	float focalDistance = 125.0f;
	float apertureSize = 1.0f;

	RTCRayHit rayIntersectScene(RTCRay ray);

	void getCosWeightedSample(Vector3 n, Vector3 &omega_i, float &pdf);
	Vector3 rotateVector(Vector3 v, Vector3 n);
};

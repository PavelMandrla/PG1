#pragma once

#include "simpleguidx11.h"
#include "corecrt_math_defines.h"
#include "vector3.h"
#include "structs.h"
#include "surface.h"
#include "camera.h"

#include <string>


class PathTracer : public SimpleGuiDX11 {
public:
	PathTracer(const int width, const int height, 
		const float fov_y, const Vector3 view_from, const Vector3 view_at,
		const char * config = "threads=0,verbose=3");
	~PathTracer();

	int InitDeviceAndScene(const char * config);

	int ReleaseDeviceAndScene();

	void LoadScene(const std::string file_name);

	Color4f get_pixel(const int x, const int y, const float t = 0.0f) override;

	int Ui();
private:
	Camera camera_;
	RTCDevice device_;
	RTCScene scene_;

	std::vector<Surface *> surfaces_;
	std::vector<Material *> materials_;

	float focalDistance = 125.0f;
	float apertureSize = 1.0f;

	RTCRay getRay(Vector3 origin, Vector3 direction);
	RTCRayHit rayIntersectScene(RTCRay ray);

	float getReflectionCoefficient(Vector3 d, Vector3 n, float ior1, float ior2);

	void getCosWeightedSample(Vector3 n, Vector3 &omega_i, float &pdf);
	Vector3 rotateVector(Vector3 v, Vector3 n);

	Vector3 getRefractedVector(Vector3 d, Vector3 n, float ior1, float ior2);
	Vector3 getReflectedVector(Vector3 d, Vector3 n);

	Color4f getReflectedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level);
	Color4f getRefractedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level);

	Color4f trace(RTCRay ray, int level, float rayIOR = 1);
};


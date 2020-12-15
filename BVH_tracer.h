#pragma once
#include "Tracer.h"
#include "BVH.h"

class BVH_tracer : public Tracer {
public:
	BVH_tracer(const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config = "threads=0,verbose=3");

	void LoadScene(const std::string file_name);

	Color4f get_pixel(const int x, const int y, const float t = 0.0f) override;
public:
	std::shared_ptr<BVH> bvh;

	virtual Color4f trace(RTCRay ray, int level, float rayIOR = 1) override;
};


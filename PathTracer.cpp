#include "stdafx.h"
#include "PathTracer.h"
#include "tutorials.h"
#include "objloader.h"
#include <omp.h>
#include <array>

PathTracer::PathTracer(const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config)
	: Tracer (width, height, fov_y, view_from, view_at, config) {
}

void PathTracer::getCosWeightedSample(Vector3 n, Vector3 & omega_i, float & pdf) {
	int threadId = omp_get_thread_num();
	float xi_1 = this->rngs[threadId].getRandNum(0, 0.99999f);
	float xi_2 = this->rngs[threadId].getRandNum(0, 0.99999f);

	Vector3 dir{
		cos(float(2 * M_PI * xi_1)) * sqrt(1.0f - xi_2),
		sin(float(2 * M_PI * xi_1)) * sqrt(1.0f - xi_2),
		sqrt(xi_2)
	};
	dir.Normalize();

	pdf = (Vector3{ 0,0,1 }.DotProduct(dir)) / M_PI;
	omega_i = this->rotateVector(dir, n);
}

Vector3 PathTracer::rotateVector(Vector3 v, Vector3 n) {
	Vector3 o1 = (abs(n.x) > abs(n.z)) ? Vector3(-n.y, n.x, 0.0f) : Vector3(0.0f, -n.z, n.y);
	Vector3 o2 = o1.CrossProduct(n);
	return Matrix3x3{ o1, o2, n } *v;
}

Color4f PathTracer::get_pixel(const int x, const int y, const float t) {
	const int multisampling_width = 20;
	const int multisamplingTotal = multisampling_width * multisampling_width;
	std::array<std::array<Color4f, multisampling_width>, multisampling_width> result_colors;

	#pragma omp parallel for num_threads(this->threadCount)
	for (int fieldX = 0; fieldX < multisampling_width; fieldX++) {
		int tid = omp_get_thread_num();
		float msX = fieldX * (1.0f / multisampling_width);
		for (int fieldY = 0; fieldY < multisampling_width; fieldY++) {
			float msY = fieldY * (1.0f / multisampling_width);
			
			float rand1 = this->rngs[tid].getRandNum(-0.5f / multisampling_width, 0.5f / multisampling_width);
			float rand2 = this->rngs[tid].getRandNum(-0.5f / multisampling_width, 0.5f / multisampling_width);

			RTCRay primaryRay = camera_.GenerateRay(x + msX + rand1, y + msY + rand2, this->focalDistance, this->apertureSize, this->rngs[tid]);
			result_colors[fieldX][fieldY] = trace(primaryRay, 0);
		}
	}

	Color4f tmpMultisamplingColor{ 0.0f, 0.0f, 0.0f, 1.0f };
	for (int fieldX = 0; fieldX < multisampling_width; fieldX++) {
		for (int fieldY = 0; fieldY < multisampling_width; fieldY++) {
			tmpMultisamplingColor.r += result_colors[fieldX][fieldY].r;
			tmpMultisamplingColor.g += result_colors[fieldX][fieldY].g;
			tmpMultisamplingColor.b += result_colors[fieldX][fieldY].b;
		}
	}
	return Color4f{ tmpMultisamplingColor.r / multisamplingTotal, tmpMultisamplingColor.g / multisamplingTotal, tmpMultisamplingColor.b / multisamplingTotal, 1.0f };
}

Color4f PathTracer::trace(RTCRay ray, int level, float rayIOR) {
	if (level > 20) {
		return Color4f{ 0, 0, 0, 0 };
	}
	RTCRayHit ray_hit = this->rayIntersectScene(ray);
	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
		RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
		Material * material = (Material *)(rtcGetGeometryUserData(geometry));
		float rho = material->reflectivity > 0.95 ? 0.95 : material->reflectivity;

		if (material->isEmittingLight()) return material->getEmittedLight();	// RETURN EMITTED LIGHT, IF EMITS LIGHT

		int tid = omp_get_thread_num();
		if (rho < this->rngs[tid].getRandNum(0.0f, 1.0f)) return Color4f{ 0, 0, 0, 1 };

		Normal3f normal;
		rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);	// get interpolated normal

		Vector3 n{ normal.x, normal.y, normal.z };
		Vector3 d{ ray.dir_x, ray.dir_y, ray.dir_z };
		n.Normalize();
		d.Normalize();
		if (n.DotProduct(d) > 0) n *= -1;	//	Kontrola orientace normaly

		Vector3 org(ray.org_x, ray.org_y, ray.org_z);
		Vector3 hitPoint = org + d * ray_hit.ray.tfar;
		Vector3 omega_i;
		float pdf;

		switch (material->type) {
		case 2:	//MIRROR
		{
			return this->getReflectedLight(hitPoint, d, n, rayIOR, 0, level);
		}
		case 3:	// MATT MATERIAL
		{
			this->getCosWeightedSample(n, omega_i, pdf);
			RTCRay secondaryRay = getRay(hitPoint, omega_i);
			Color4f L_i = trace(secondaryRay, level + 1, rayIOR);
			float f_r = float(material->reflectivity / M_PI);
			//return L_i * f_r * (omega_i.DotProduct(n))
			float tmp = f_r * (omega_i.DotProduct(n)) / (pdf * rho);

			return Color4f{ 
				L_i.r * material->diffuse.x * tmp,
				L_i.g * material->diffuse.y * tmp,
				L_i.b * material->diffuse.z * tmp,
				1.0f
			};
		}

		case 4:	// DIELECTRIC MATERIAL
		{
			float matIOR = rayIOR == material->ior ? 1 : material->ior;
			auto reflectedLight = this->getReflectedLight(hitPoint, d, n, rayIOR, matIOR, level);
			auto refractedLight = this->getRefractedLight(hitPoint, d, n, rayIOR, matIOR, level);

			float attenuationRed, attenuationGreen, attenuationBlue;
			if (rayIOR == material->ior) {	// inside material
				attenuationRed = exp(-(1 - material->diffuse.x)*ray.tfar);
				attenuationGreen = exp(-(1 - material->diffuse.y)*ray.tfar);
				attenuationBlue = exp(-(1 - material->diffuse.z)*ray.tfar);
			} else {						// in air
				attenuationRed = 1;
				attenuationGreen = 1;
				attenuationBlue = 1;
			}

			float r = this->getReflectionCoefficient(d, n, rayIOR, matIOR);
			Color3f attenuation{ attenuationRed, attenuationGreen, attenuationBlue };
			reflectedLight.a = r;
			refractedLight.a = 1.0f - r;

			return (reflectedLight + refractedLight) * attenuation;
		}
		}

	}
	//return this->sphericalMap->getBackgroundColor(ray.dir_x, ray.dir_y, ray.dir_z);
	return Color4f{ 0,0,0,1 };
}
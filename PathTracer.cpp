#include "stdafx.h"
#include "PathTracer.h"
#include "tutorials.h"
#include "objloader.h"
#include <omp.h>
#include <array>

PathTracer::PathTracer(const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config) : SimpleGuiDX11(width, height) {
	InitDeviceAndScene(config);
	camera_ = Camera(width, height, fov_y, view_from, view_at);
}

PathTracer::~PathTracer() {
	ReleaseDeviceAndScene();
}

int PathTracer::InitDeviceAndScene(const char * config) {
	device_ = rtcNewDevice(config);
	error_handler(nullptr, rtcGetDeviceError(device_), "Unable to create a new device.\n");
	rtcSetDeviceErrorFunction(device_, error_handler, nullptr);

	ssize_t triangle_supported = rtcGetDeviceProperty(device_, RTC_DEVICE_PROPERTY_TRIANGLE_GEOMETRY_SUPPORTED);

	// create a new scene bound to the specified device
	scene_ = rtcNewScene(device_);

	return S_OK;
}

int PathTracer::ReleaseDeviceAndScene() {
	rtcReleaseScene(scene_);
	rtcReleaseDevice(device_);

	return S_OK;
}

void PathTracer::LoadScene(const std::string file_name) {
	const int no_surfaces = LoadOBJ(file_name.c_str(), surfaces_, materials_);

	// surfaces loop
	for (auto surface : surfaces_) {
		RTCGeometry mesh = rtcNewGeometry(device_, RTC_GEOMETRY_TYPE_TRIANGLE);

		Vertex3f * vertices = (Vertex3f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof(Vertex3f), 3 * surface->no_triangles());

		Triangle3ui * triangles = (Triangle3ui *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof(Triangle3ui), surface->no_triangles());

		rtcSetGeometryUserData(mesh, (void*)(surface->get_material()));

		rtcSetGeometryVertexAttributeCount(mesh, 2);

		Normal3f * normals = (Normal3f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof(Normal3f), 3 * surface->no_triangles());

		Coord2f * tex_coords = (Coord2f *)rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof(Coord2f), 3 * surface->no_triangles());


		// triangles loop
		for (int i = 0, k = 0; i < surface->no_triangles(); ++i) {
			Triangle & triangle = surface->get_triangle(i);
			triangle.material = surface->get_material();

			// vertices loop
			for (int j = 0; j < 3; ++j, ++k) {
				const Vertex & vertex = triangle.vertex(j);

				vertices[k].x = vertex.position.x;
				vertices[k].y = vertex.position.y;
				vertices[k].z = vertex.position.z;

				normals[k].x = vertex.normal.x;
				normals[k].y = vertex.normal.y;
				normals[k].z = vertex.normal.z;

				tex_coords[k].u = vertex.texture_coords[0].u;
				tex_coords[k].v = vertex.texture_coords[0].v;
			} // end of vertices loop

			triangles[i].v0 = k - 3;
			triangles[i].v1 = k - 2;
			triangles[i].v2 = k - 1;
		} // end of triangles loop

		rtcCommitGeometry(mesh);
		unsigned int geom_id = rtcAttachGeometry(scene_, mesh);
		rtcReleaseGeometry(mesh);
	} // end of surfaces loop

	rtcCommitScene(scene_);
}

int PathTracer::Ui() {
	static float f = 0.0f;
	static int counter = 0;

	// Use a Begin/End pair to created a named window
	ImGui::Begin("Ray Tracer Params");

	ImGui::Text("Surfaces = %d", surfaces_.size());
	ImGui::Text("Materials = %d", materials_.size());
	ImGui::Separator();
	ImGui::Checkbox("Vsync", &vsync_);

	if (ImGui::Button("Button"))
		counter++;
	ImGui::SameLine();
	ImGui::Text("counter = %d", counter);

	ImGui::SliderFloat("aperture", &this->apertureSize, 0.0f, 5.0f); // Edit 1 float using a slider from 0.0f to 1.0f   
	ImGui::SliderFloat("focal point", &this->focalDistance, 50.0f, 300.0f); // Edit 1 float using a slider from 0.0f to 1.0f   

	ImGui::Text("Application average %.3f ms/frame (%.1f FPS)", 1000.0f / ImGui::GetIO().Framerate, ImGui::GetIO().Framerate);
	ImGui::End();

	return 0;
}

RTCRay PathTracer::getRay(Vector3 origin, Vector3 direction) {
	RTCRay ray;
	ray.org_x = origin.x;
	ray.org_y = origin.y;
	ray.org_z = origin.z;

	ray.dir_x = direction.x;
	ray.dir_y = direction.y;
	ray.dir_z = direction.z;

	ray.tnear = 0.0001f;
	ray.tfar = FLT_MAX;
	ray.time = 0.0f;

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

float PathTracer::getReflectionCoefficient(Vector3 d, Vector3 n, float ior1, float ior2) {
	float r0 = pow((ior1 - ior2) / (ior1 + ior2), 2);
	float theta_t = acos(n.DotProduct(d));
	float theta_i = asin((ior1 * sin(theta_t)) / ior2);
	float theta = ior1 <= ior2 ? theta_i : theta_t;
	if (theta > M_PI_2) theta = float(M_PI - theta);
	return r0 + (1 - r0)*pow(1 - cos(theta), 5);
}

RTCRayHit PathTracer::rayIntersectScene(RTCRay ray) {
	RTCHit hit;
	hit.geomID = RTC_INVALID_GEOMETRY_ID;
	hit.primID = RTC_INVALID_GEOMETRY_ID;
	hit.Ng_x = 0.0f; // geometry normal
	hit.Ng_y = 0.0f;
	hit.Ng_z = 0.0f;

	// merge ray and hit structures
	RTCRayHit ray_hit;
	ray_hit.ray = ray;
	ray_hit.hit = hit;

	// intersect ray with the scene
	RTCIntersectContext context;
	rtcInitIntersectContext(&context);
	rtcIntersect1(scene_, &context, &ray_hit);

	return ray_hit;
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


Vector3 PathTracer::getRefractedVector(Vector3 d, Vector3 n, float ior1, float ior2) {
	float iorRatio = ior1 / ior2;
	float phi1 = acos(Vector3{ d.x*-1, d.y*-1, d.z*-1 }.DotProduct(n));
	float phi2 = asin(sqrt(1 - pow(cos(phi1), 2)));

	return iorRatio * d + (iorRatio*cos(phi1) - cos(phi2)) * n;
}

Vector3 PathTracer::getReflectedVector(Vector3 d, Vector3 n) {
	d.Normalize();
	n.Normalize();
	Vector3 vec = d - 2 * (d.DotProduct(n))*n;
	return vec;
}


Color4f PathTracer::getReflectedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level) {
	if (level > 10) return Color4f{ 0, 0, 0, 1.0f };
	Vector3 reflectedVector = getReflectedVector(d, n);
	RTCRay reflectedRay = getRay(hitPoint, reflectedVector);
	return trace(reflectedRay, level + 1, envIOR);
}

Color4f PathTracer::getRefractedLight(Vector3 hitPoint, Vector3 d, Vector3 n, float envIOR, float matIOR, int level) {
	Vector3 refractedVector = getRefractedVector(d, n, envIOR, matIOR);
	RTCRay refractedRay = getRay(hitPoint, refractedVector);
	return trace(refractedRay, level + 1, matIOR);
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
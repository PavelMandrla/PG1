#include "stdafx.h"
#include "LightSource.h"


LightSource::LightSource(Vector3 position, Vector3 ambient, Vector3 diffuse, Vector3 specular) {
	this->postion = position;
	this->ambient = ambient;
	this->diffuse = diffuse;
	this->specular = specular;
}


LightSource::~LightSource() {
}

RTCRay LightSource::GenerateRay(const float x_pos, const float y_pos, const float z_pos) {
	RTCRay ray = RTCRay();

	ray.org_x = this->postion.x;
	ray.org_y = this->postion.y;
	ray.org_z = this->postion.z;

	ray.dir_x = x_pos;
	ray.dir_y = y_pos;
	ray.dir_z = z_pos;

	Vector3 to(x_pos, y_pos, z_pos);
	ray.tnear = 0.01;
	ray.tfar = (this->postion - to).L2Norm();
	ray.time = 0.0f;

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

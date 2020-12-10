#include "stdafx.h"
#include "camera.h"
#include <cmath>

Camera::Camera( const int width, const int height, const float fov_y,
	const Vector3 view_from, const Vector3 view_at )
{
	width_ = width;
	height_ = height;
	fov_y_ = fov_y;

	view_from_ = view_from;
	view_at_ = view_at;

	// TODO compute focal lenght based on the vertical field of view and the camera resolution
	f_y_ = height_ / (2 * tan(fov_y_ / 2));
	// TODO build M_c_w_ matrix	
	Vector3 z_c = view_from_ - view_at_;
	Vector3 x_c = up_.CrossProduct(z_c);
	Vector3 y_c = z_c.CrossProduct(x_c);

	x_c.Normalize();
	y_c.Normalize();
	z_c.Normalize();

	M_c_w_ = Matrix3x3(x_c, y_c, z_c);

	std::default_random_engine(std::random_device{}());
}

RTCRay Camera::GenerateRay( const float x_i, const float y_i ) const
{
	RTCRay ray = RTCRay();

	// TODO fill in ray structure and compute ray direction
	ray.org_x = view_from_.x;
	ray.org_y = view_from_.y;
	ray.org_z = view_from_.z;
	
	Vector3 d_c(x_i - width_ / 2, height_ / 2 - y_i, -f_y_);
	d_c.Normalize();
	Vector3 d_c_ws = M_c_w_ * d_c;

	ray.dir_x = d_c_ws.x;
	ray.dir_y = d_c_ws.y;
	ray.dir_z = d_c_ws.z;

	ray.tnear = FLT_MIN;
	ray.tfar = FLT_MAX;
	ray.time = 0.0f;

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

RTCRay Camera::GenerateRay(const float x_i, const float y_i, const float focalDistance, const float apertureSize, RNG &rng) const {

	Vector3 d_c(x_i - width_ / 2, height_ / 2 - y_i, -f_y_);
	d_c.Normalize();
	Vector3 d_c_ws = M_c_w_ * d_c;
	d_c_ws.Normalize();
	Vector3 focalPoint = Vector3 { view_from_.x, view_from_.y , view_from_.z } + Vector3{ d_c_ws } * focalDistance;

	RTCRay ray = RTCRay();
	
	Vector3 originShift = M_c_w_ * Vector3 {
		rng.getRandNum(-apertureSize / 2.0f, apertureSize / 2.0f),
		rng.getRandNum(-apertureSize / 2.0f, apertureSize / 2.0f),
		0
	};

	ray.org_x = view_from_.x + originShift.x;
	ray.org_y = view_from_.y + originShift.y;
	ray.org_z = view_from_.z + originShift.z;

	Vector3 nDir = focalPoint - Vector3{ ray.org_x, ray.org_y , ray.org_z };
	nDir.Normalize();

	ray.dir_x = nDir.x;
	ray.dir_y = nDir.y;
	ray.dir_z = nDir.z;
	
	ray.tnear = FLT_MIN;
	ray.tfar = FLT_MAX;
	ray.time = 0.0f;

	ray.mask = 0; // can be used to mask out some geometries for some rays
	ray.id = 0; // identify a ray inside a callback function
	ray.flags = 0; // reserved

	return ray;
}

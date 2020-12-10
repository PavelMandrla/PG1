#include "stdafx.h"
#include "SphericalMap.h"
#include "corecrt_math_defines.h"
#include "vector3.h"

#include <iostream>

SphericalMap::SphericalMap() {
}

SphericalMap::SphericalMap(const std::string & file_name) {
	this->texture_ = std::shared_ptr<Texture>(new Texture(file_name.c_str()));
	//this->texture_ = std::make_unique<Texture>(new Texture(file_name.c_str()));
}

Color3f SphericalMap::texel(const float x, const float y, const float z) const {
	Vector3 vec{ x, y, z };
	vec.Normalize();
	float u = 0.5 + atan2(vec.x, vec.z) / (2 * M_PI);
	float v = 0.5 - asin(vec.y) / M_PI;
	return this->texture_->get_texel(u, v);
}

Color4f SphericalMap::getBackgroundColor(const float x, const float y, const float z) const {
	auto backgroundColor = this->texel(x, z, y);
	return Color4f{ backgroundColor.r, backgroundColor.g, backgroundColor.b, 1.0f };
}

SphericalMap::~SphericalMap() {
}

#pragma once

#include <string>
#include "texture.h"

class SphericalMap {
public:
	SphericalMap();

	SphericalMap( const std::string & file_name );

	Color3f texel(const float x, const float y, const float z) const;

	Color4f getBackgroundColor(const float x, const float y, const float z) const;

	~SphericalMap();
private:
	std::shared_ptr<Texture> texture_;

};


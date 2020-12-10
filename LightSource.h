#pragma once

#include "vector3.h"

class LightSource {
public:
	Vector3 postion;

	Vector3 ambient; /*!< RGB barva prostředí \f$\left<0, 1\right>^3\f$. */
	Vector3 diffuse; /*!< RGB barva rozptylu \f$\left<0, 1\right>^3\f$. */
	Vector3 specular; /*!< RGB barva odrazu \f$\left<0, 1\right>^3\f$. */

	LightSource(Vector3 position, Vector3 ambient, Vector3 diffuse, Vector3 specular);
	~LightSource();

	RTCRay GenerateRay(const float x_i, const float y_i, const float z_i);

};


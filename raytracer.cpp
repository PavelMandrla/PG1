#include "stdafx.h"
#include "raytracer.h"
#include "objloader.h"
#include "tutorials.h"
#include <cmath>
#include <iostream>

Raytracer::Raytracer( const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config )
	: Tracer(width, height, fov_y, view_from, view_at, config) {
}

void Raytracer::LoadScene( const std::string file_name ) {
	const int no_surfaces = LoadOBJ( file_name.c_str(), surfaces_, materials_ );
	//std::vector<std::shared_ptr<Triangle>> myTriangles;	//BVH

	// surfaces loop
	for ( auto surface : surfaces_ ) {
		RTCGeometry mesh = rtcNewGeometry( device_, RTC_GEOMETRY_TYPE_TRIANGLE );

		Vertex3f * vertices = ( Vertex3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX, 0, RTC_FORMAT_FLOAT3,
			sizeof( Vertex3f ), 3 * surface->no_triangles() );

		Triangle3ui * triangles = ( Triangle3ui * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_INDEX, 0, RTC_FORMAT_UINT3,
			sizeof( Triangle3ui ), surface->no_triangles() );

		rtcSetGeometryUserData( mesh, ( void* )( surface->get_material() ) );

		rtcSetGeometryVertexAttributeCount( mesh, 2 );

		Normal3f * normals = ( Normal3f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, RTC_FORMAT_FLOAT3,
			sizeof( Normal3f ), 3 * surface->no_triangles() );

		Coord2f * tex_coords = ( Coord2f * )rtcSetNewGeometryBuffer(
			mesh, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, RTC_FORMAT_FLOAT2,
			sizeof( Coord2f ), 3 * surface->no_triangles() );		

		
		// triangles loop
		for ( int i = 0, k = 0; i < surface->no_triangles(); ++i )
		{
			Triangle & triangle = surface->get_triangle( i );
			triangle.material = surface->get_material();
			//myTriangles.push_back(std::make_shared<Triangle>(triangle));	//BVH

			// vertices loop
			for ( int j = 0; j < 3; ++j, ++k )
			{
				const Vertex & vertex = triangle.vertex( j );

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

		rtcCommitGeometry( mesh );
		unsigned int geom_id = rtcAttachGeometry( scene_, mesh );
		rtcReleaseGeometry( mesh );
	} // end of surfaces loop
	/*
	this->bvh = std::make_shared<BVH>(myTriangles);		//BVH
	this->bvh->initialize();								//BVH
	*/

	//this->lightSources.push_back(LightSource(Vector3(-200, 600, -30), Vector3(0.8f, 0.8f, 0.8f), Vector3(1, 1, 1), Vector3(1, 1, 1)));
	this->lightSources.push_back(LightSource(Vector3(100,0, 300), Vector3(0.8, 0.8f, 0.8f), Vector3(1, 1, 1), Vector3(1, 1, 1)));
	//this->lightSources.push_back(LightSource(Vector3(175, -140, 130), Vector3(0.8, 0.8f, 0.8f), Vector3(1, 1, 1), Vector3(1, 1, 1)));
	//this->lightSources.push_back(LightSource(Vector3(175, -140, 130), Vector3(0.8f, 0.8f, 0.8f), Vector3(1, 1, 1), Vector3(1, 1, 1)));
	this->ambient = new LightSource(Vector3(0, 0, 0), Vector3(0.1f, 0.1f, 0.1f), Vector3(0, 0, 0), Vector3(0.0f, 0.0f, 0.0f));
	this->sphericalMap = std::shared_ptr<SphericalMap>(new SphericalMap("../../../data/museumplein.jpg"));

	rtcCommitScene( scene_ );
}

Color4f Raytracer::trace(RTCRay ray, int level, float rayIOR) {
	if (level > 6) return Color4f{ 0,0,0,1 };
	RTCRayHit ray_hit = this->rayIntersectScene(ray);

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID) { // we hit something
		RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
		Normal3f normal;
		// get interpolated normal
		rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 0, &normal.x, 3);

		Vector3 n { normal.x, normal.y, normal.z };		//	Kontrola orientace normaly
		Vector3 d { ray.dir_x, ray.dir_y, ray.dir_z };
		if (n.DotProduct(d) > 0) n *= -1;
		n.Normalize();
		d.Normalize();

		Material * material = (Material *)(rtcGetGeometryUserData(geometry));
		
		Vector3 org(ray.org_x, ray.org_y, ray.org_z);
		Vector3 hitPoint = org + d * ray_hit.ray.tfar;

		switch (material->type) {
			case 1: { // LAMBERT SHADER
				Color4f result{ 0,0,0,1 };
				auto C = Color4f {
					material->diffuse.x,
					material->diffuse.y,
					material->diffuse.z,
					1.0f
				};

				for (auto lightSource : this->lightSources) {
					if (!this->isBlocked(lightSource, hitPoint)) {
						Vector3 L = lightSource.postion - hitPoint;
						L.Normalize();
						auto LN = L.DotProduct(n);

						result += Color4f{
							LN * C.r * lightSource.diffuse.x,
							LN * C.g * lightSource.diffuse.y,
							LN * C.b * lightSource.diffuse.z,
							1
						};
					}
				}
				return result;
			}

			case 2: { // NORMAL SHADER
				float r_n = (n.x + 1) / 2.0f;
				float g_n = (n.y + 1) / 2.0f;
				float b_n = (n.z + 1) / 2.0f;
				return Color4f{ b_n, g_n, r_n, 1.0f };
			}

			case 3: { //PHONG SHADER		
				float matIOR = rayIOR == material->ior ? 1 : material->ior;
				auto phong = getPhongIllumination(hitPoint, d, n, material, geometry, ray_hit);
				auto reflectedLight = this->getReflectedLight(hitPoint, d, n, rayIOR, matIOR, level);
				float r = this->getReflectionCoefficient(d, n, rayIOR, matIOR);
				reflectedLight.a = r;
				phong.a = 1.0f - r;
				return phong + reflectedLight;
			}

			case 4: { // DIELECTRIC SHADER					
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
	return this->sphericalMap->getBackgroundColor(ray.dir_x, ray.dir_y, ray.dir_z);
}

Color4f Raytracer::getPhongIllumination(Vector3 hitPoint, Vector3 d, Vector3 n, Material *material, RTCGeometry geometry, RTCRayHit ray_hit) {
	float ligtVal_R = ambient->ambient.x * material->ambient.x;
	float ligtVal_G = ambient->ambient.y * material->ambient.y;
	float ligtVal_B = ambient->ambient.z * material->ambient.z;

	float diffuseR, diffuseG, diffuseB;
	Texture * diffuseTexture = material->get_texture(Material::kDiffuseMapSlot);
	if (diffuseTexture == nullptr) {
		diffuseR = material->diffuse.x;
		diffuseG = material->diffuse.y;
		diffuseB = material->diffuse.z;
	} else {
		Coord2f tex_coord;
		rtcInterpolate0(geometry, ray_hit.hit.primID, ray_hit.hit.u, ray_hit.hit.v, RTC_BUFFER_TYPE_VERTEX_ATTRIBUTE, 1, &tex_coord.u, 2);
		auto texel = diffuseTexture->get_texel(tex_coord.u, 1.0f - tex_coord.v);
		diffuseR = texel.r;
		diffuseG = texel.g;
		diffuseB = texel.b;
	}
	
	
	for (auto lightSrc : this->lightSources) {
		if (!isBlocked(lightSrc, hitPoint)) {
			Vector3 l = lightSrc.postion - hitPoint;
			l.Normalize();
			Vector3 lr = 2 * (n.DotProduct(l)) * n - l;

			ligtVal_R += lightSrc.diffuse.x * diffuseR * (n.DotProduct(l));
			ligtVal_G += lightSrc.diffuse.y * diffuseG * (n.DotProduct(l));
			ligtVal_B += lightSrc.diffuse.z * diffuseB * (n.DotProduct(l));

			ligtVal_R += lightSrc.specular.x * material->specular.x * pow((-1 * d).DotProduct(lr), material->shininess);
			ligtVal_G += lightSrc.specular.y * material->specular.y * pow((-1 * d).DotProduct(lr), material->shininess);
			ligtVal_B += lightSrc.specular.z * material->specular.z * pow((-1 * d).DotProduct(lr), material->shininess);
		}
	}
	return Color4f{ ligtVal_R, ligtVal_G, ligtVal_B, 1.0 };

}

bool Raytracer::isBlocked(LightSource lightSrc, Vector3 hitPoint) {
	RTCRay ray = lightSrc.GenerateRay(hitPoint.x, hitPoint.y, hitPoint.z);

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

	if (ray_hit.hit.geomID != RTC_INVALID_GEOMETRY_ID) {
		RTCGeometry geometry = rtcGetGeometry(scene_, ray_hit.hit.geomID);
		Material * material = (Material *)(rtcGetGeometryUserData(geometry));
		if (material->type == 4) {
			Vector3 org{ ray.org_x, ray.org_y, ray.org_z };
			Vector3 d { ray.dir_x, ray.dir_y, ray.dir_z };
			Vector3 new_hitPoint = org + d * ray_hit.ray.tfar;
			return isBlocked(lightSrc, new_hitPoint);
		}
		return true;
	}
	return false;
}

Color4f Raytracer::get_pixel( const int x, const int y, const float t ) {
	
	const int multisampling_width = 3;
	const int multisamplingTotal = multisampling_width * multisampling_width;

	std::array<std::array<Color4f, multisampling_width>, multisampling_width> result_colors;

	if (multisampling_width == 1) {
		RTCRay primaryRay = camera_.GenerateRay(x, y);
		return trace(primaryRay, 0);
	}

	#pragma omp parallel for num_threads(this->threadCount)
	for (int fieldX = 0; fieldX < multisampling_width; fieldX++) {
		
		float msX = fieldX * (1.0f / multisampling_width);
		for (int fieldY = 0; fieldY < multisampling_width; fieldY++) {
			float msY = fieldY * (1.0f / multisampling_width);
			int tid = omp_get_thread_num();

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
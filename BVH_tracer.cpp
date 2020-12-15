#include "stdafx.h"
#include "BVH_tracer.h"
#include "objloader.h"

BVH_tracer::BVH_tracer(const int width, const int height, const float fov_y, const Vector3 view_from, const Vector3 view_at, const char * config)
	: Tracer(width, height, fov_y, view_from, view_at, config) {
}

void BVH_tracer::LoadScene(const std::string file_name) {
	const int no_surfaces = LoadOBJ(file_name.c_str(), surfaces_, materials_);
	std::vector<Triangle*> myTriangles;	//BVH

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
			myTriangles.push_back(&triangle);	//BVH

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
	
	this->bvh = std::make_shared<BVH>(myTriangles);		//BVH
	this->bvh->initialize();								//BVH

	rtcCommitScene(scene_);
}

Color4f BVH_tracer::get_pixel(const int x, const int y, const float t) {
	RTCRay primaryRay = camera_.GenerateRay(x, y);
	return trace(primaryRay, 0);
}

Color4f BVH_tracer::trace(RTCRay ray, int level, float rayIOR) {
	std::shared_ptr<Triangle> triangle = this->bvh->traverse(ray);	//BVH
	if (triangle != nullptr) {
		return Color4f{
			triangle->material->diffuse.x,
			triangle->material->diffuse.y,
			triangle->material->diffuse.z,
			1.0f
		};
	}
	return Color4f{ 0.1f, 0.1f ,0.1f ,1.0f };
}

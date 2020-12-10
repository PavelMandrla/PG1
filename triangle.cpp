#include "stdafx.h"
#include "triangle.h"

Triangle::Triangle( const Vertex & v0, const Vertex & v1, const Vertex & v2, Surface * surface )
{
	vertices_[0] = v0;
	vertices_[1] = v1;
	vertices_[2] = v2;	

	// ukazatel na surface schováme (!pokud se tam vejde!) do paddingu prvního vertexu
	*reinterpret_cast<Surface **>( &vertices_[0].pad ) = surface;	
}

Vertex Triangle::vertex( const int i ) {
	return vertices_[i];
}

Surface * Triangle::surface()
{	
	return *reinterpret_cast<Surface **>( vertices_[0].pad ); // FIX: chybí verze pro 64bit
}

Vector3 Triangle::getCenter() {
	Vector3 center;
	center.x = (vertices_[0].position.x + vertices_[1].position.x + vertices_[2].position.x) / 3;
	center.y = (vertices_[0].position.y + vertices_[1].position.y + vertices_[2].position.y) / 3;
	center.z = (vertices_[0].position.z + vertices_[1].position.z + vertices_[2].position.z) / 3;
	return center;
}

Vector3 Triangle::getNormal() {
	Vector3 normal = (vertices_[1].position - vertices_[0].position).CrossProduct(vertices_[2].position - vertices_[0].position);
	normal.Normalize();
	return normal;
}

std::shared_ptr<Vector3> Triangle::intersect(RTCRay ray) { //PODLE -> https://www.scratchapixel.com/lessons/3d-basic-rendering/ray-tracing-rendering-a-triangle/ray-triangle-intersection-geometric-solution
	Vector3 rayDir { ray.dir_x, ray.dir_y, ray.dir_z };
	Vector3 rayOrig { ray.org_x, ray.org_y, ray.org_z };
	Vector3 n = this->getNormal();
	
	if (n.DotProduct(rayDir) == 0) return nullptr;		//triangle and ray are paralel -> no intersection

	float d = n.DotProduct(vertices_[0].position);
	float t = (n.DotProduct(rayOrig) + d) / n.DotProduct(rayDir);
	if (t < 0) return nullptr;

	Vector3 p = rayOrig + t * rayDir;	//nalezení průsečíku s rovinou

	//kontrola, zda se průsečík s rovinou definovanou trojúhelníkem nachází uvnitř
	Vector3 c;

	Vector3 edge0 = vertices_[1].position - vertices_[0].position;
	Vector3 vp0 = p - vertices_[0].position;
	c = edge0.CrossProduct(vp0);
	if (n.DotProduct(c) < 0) return nullptr;

	Vector3 edge1 = vertices_[2].position - vertices_[1].position;
	Vector3 vp1 = p - vertices_[1].position;
	c = edge1.CrossProduct(vp1);
	if (n.DotProduct(c) < 0) return nullptr;

	Vector3 edge2 = vertices_[0].position - vertices_[2].position;
	Vector3 vp2 = p - vertices_[2].position;
	c = edge2.CrossProduct(vp2);
	if (n.DotProduct(c) < 0) return nullptr;

	return std::make_shared<Vector3>(p);
}


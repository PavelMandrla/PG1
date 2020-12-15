#include "stdafx.h"
#include "BVH.h"

AABB::AABB(Node* node) {
	auto bvh = node->bvh.lock();
	std::vector<float> x_points, y_points, z_points;

	for (int i = node->from; i <= node->to; i++) {
		for (int j = 0; j < 3; j++) {
			x_points.push_back(bvh->triangles[i]->vertex(j).position.x);
			y_points.push_back(bvh->triangles[i]->vertex(j).position.y);
			z_points.push_back(bvh->triangles[i]->vertex(j).position.z);
		}
	}

	this->bounds[0] = Vector3{
		*std::min_element(x_points.begin(), x_points.end()),
		*std::min_element(y_points.begin(), y_points.end()),
		*std::min_element(z_points.begin(), z_points.end())
	};
	this->bounds[1] = Vector3{
		*std::max_element(x_points.begin(), x_points.end()),
		*std::max_element(y_points.begin(), y_points.end()),
		*std::max_element(z_points.begin(), z_points.end())
	};
}

bool AABB::isIntersected(RTCRay ray) {
	Vector3 d = {
		ray.dir_x,
		ray.dir_y,
		ray.dir_z
	};
	d.Normalize();

	std::vector<float> t_0 {
		(this->bounds[0].x - ray.org_x) / d.x,
		(this->bounds[0].y - ray.org_y) / d.y,
		(this->bounds[0].z - ray.org_z) / d.z
	};
	std::vector<float> t_1 {
		(this->bounds[1].x - ray.org_x) / d.x,
		(this->bounds[1].y - ray.org_y) / d.y,
		(this->bounds[1].z - ray.org_z) / d.z
	};

	if (t_0[0] > t_1[0]) std::swap(t_0[0], t_1[0]);
	if (t_0[1] > t_1[1]) std::swap(t_0[1], t_1[1]);
	if (t_0[2] > t_1[2]) std::swap(t_0[2], t_1[2]);

	float t0 = *(std::max_element)(t_0.begin(), t_0.end());
	float t1 = *(std::min_element)(t_1.begin(), t_1.end());

	if ((t0 <= t1) && (0 < t1)) {
		return true;
	}
	return false;
}

float AABB::getSurfaceArea() {
	Vector3 d = this->bounds[1] - this->bounds[0];
	return 2 * (d.x * d.y + d.x * d.z + d.y * d.z);
}

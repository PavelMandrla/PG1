#include "stdafx.h"
#include "BVH.h"

BVH::Node::Node(int from, int to, std::shared_ptr<BVH> bvh) {
	this->span[0] = from;
	this->span[1] = to;
	this->children[0] = children[1] = nullptr;
	this->bounds = nullptr;
	this->bvh = bvh;
}

std::shared_ptr<AABB> BVH::Node::getAABB() {
	if (this->bounds == nullptr) {
		std::vector<float> pointsX, pointsY, pointsZ;
		if (this->isLeaf()) {
			auto bvh_sp = bvh.lock();
			auto v0 = bvh_sp->items.at(span[0])->vertex(0);
			auto v1 = bvh_sp->items.at(span[0])->vertex(1);
			auto v2 = bvh_sp->items.at(span[0])->vertex(2);

			pointsX = { v0.position.x, v1.position.x, v2.position.x };
			pointsY = { v0.position.y, v1.position.y, v2.position.y };
			pointsZ = { v0.position.z, v1.position.z, v2.position.z };
		} else {
			auto ch0AABB = children[0]->getAABB();
			auto ch1AABB = children[1]->getAABB();

			pointsX = { ch0AABB->bounds[0].x, ch0AABB->bounds[1].x, ch1AABB->bounds[0].x, ch1AABB->bounds[1].x };
			pointsY = { ch0AABB->bounds[0].y, ch0AABB->bounds[1].y, ch1AABB->bounds[0].y, ch1AABB->bounds[1].y };
			pointsZ = { ch0AABB->bounds[0].z, ch0AABB->bounds[1].z, ch1AABB->bounds[0].z, ch1AABB->bounds[1].z };
		}
		this->bounds = std::make_shared<AABB>(
			Vector3 {
				*std::min_element(pointsX.begin(), pointsX.end()),
				*std::min_element(pointsY.begin(), pointsY.end()),
				*std::min_element(pointsZ.begin(), pointsZ.end())
			},
			Vector3 {
				*std::max_element(pointsX.begin(), pointsX.end()),
				*std::max_element(pointsY.begin(), pointsY.end()),
				*std::max_element(pointsZ.begin(), pointsZ.end())
			}
		);
		
	}
	return this->bounds;
}

TriangleAndDistance BVH::Node::intersect(RTCRay & ray) {
	if (this->isIntersected(ray)) {
		if (this->isLeaf()) {
			auto bvh_sp = bvh.lock();
			std::shared_ptr<Vector3> intersection  = bvh_sp->items.at(this->span[0])->intersect(ray);
			if (intersection == nullptr) return TriangleAndDistance{ nullptr, FLT_MAX };

			return TriangleAndDistance{
				bvh_sp->items.at(this->span[0]),
				Vector3{ intersection->x - ray.org_x, intersection->y - ray.org_y, intersection->z - ray.org_z }.L2Norm()
			};
		} else {
			TriangleAndDistance t1 = this->children[0]->intersect(ray);
			TriangleAndDistance t2 = this->children[1]->intersect(ray);

			return t1.distance < t2.distance ? t1 : t2;	//return closer triangle
		}
	}
	return TriangleAndDistance{ nullptr, FLT_MAX };
}

bool BVH::Node::isIntersected(RTCRay & ray) {
	float tx0 = (getAABB()->bounds[0].x - ray.org_x) / ray.dir_x;
	float ty0 = (getAABB()->bounds[0].y - ray.org_y) / ray.dir_y;
	float tz0 = (getAABB()->bounds[0].z - ray.org_z) / ray.dir_z;
	
	float tx1 = (getAABB()->bounds[1].x - ray.org_x) / ray.dir_x;
	float ty1 = (getAABB()->bounds[1].y - ray.org_y) / ray.dir_y;
	float tz1 = (getAABB()->bounds[1].z - ray.org_z) / ray.dir_z;

	if (tx0 > tx1) std::swap(tx0, tx1);
	if (ty0 > ty1) std::swap(ty0, ty1);
	if (tz0 > tz1) std::swap(tz0, tz1);

	std::vector<float> tmpT0 = { tx0, ty0, tz0 };
	std::vector<float> tmpT1 = { tx1, ty1, tz1 };

	float t0 = *std::max_element(tmpT0.begin(), tmpT0.end());
	float t1 = *std::min_element(tmpT1.begin(), tmpT1.end());

	if (t0 < t1) {
		return true;
	}
	
	return false;
}

bool BVH::Node::isLeaf() {
	return this->span[0] == this->span[1];
}

BVH::BVH(std::vector<std::shared_ptr<Triangle>> items) {
	this->items = items;
}

BVH::~BVH() {}


std::shared_ptr<BVH::Node> BVH::BuildTree(int from, int to, int depth) {
	//std::cout << from << " -> {" << from << ", " << to << "}" << std::endl;
	auto result = std::make_shared<Node>(from, to, shared_from_this()); //problem
	if (from == to) return result;

	std::sort(this->items.begin() + from, this->items.begin() + to, [depth](const std::shared_ptr<Triangle> &a, const std::shared_ptr<Triangle> &b) -> bool {
		switch (depth % 3) {
		case 0:
			return a->getCenter().x > b->getCenter().x;
		case 1:
			return a->getCenter().y > b->getCenter().y;
		case 2:
			return a->getCenter().z > b->getCenter().z;
		}
	});

	int pivotIndex = floor(from + (to - from) / 2);
	result->children[0] = BuildTree(from, pivotIndex, depth+1);
	result->children[1] = BuildTree(pivotIndex+1, to, depth+1);
	
	return result;
}

void BVH::initialize() {
	this->root = this->BuildTree(0, items.size() - 1, 0);
	this->root->getAABB();
}

std::shared_ptr<Triangle> BVH::traverse(RTCRay &ray) {
	TriangleAndDistance res = this->root->intersect(ray);
	ray.tfar = res.distance;
	return res.triangle;
}




#include "stdafx.h"
#include "BVH.h"

BVH::BVH(std::vector<Triangle*> triangles) {
	this->triangles = triangles;
}

void BVH::initialize() {
	this->root = std::make_shared<Node>(0, this->triangles.size() -1, shared_from_this());
	this->root->buildTree(0);
}

std::shared_ptr<Triangle> BVH::traverse(RTCRay &ray) {
	return this->root->traverse(ray);
}

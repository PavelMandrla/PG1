#include "stdafx.h"
#include "BVH.h"
#include <vector>
#include <algorithm>

#include <iostream>

Node::Node(int from, int to, std::shared_ptr<BVH> bvh) {
	this->from = from;
	this->to = to;
	this->bvh = bvh;
	this->isLeaf_ = false;
}

bool Node::isLeaf() {
	return this->isLeaf_;
}

std::shared_ptr<AABB> Node::getAABB() {
	if (this->aabb == nullptr) {
		this->aabb = std::make_shared<AABB>(this);
	}
	return this->aabb;
}


float Node::calculateSAH(Division div) {
	float SA = this->getAABB()->getSurfaceArea();
	const float t_traversal = 1.0f;
	const float t_intersect = 2.0f;

	return t_traversal + div.left->getSAHNodeVal(SA) + div.right->getSAHNodeVal(SA);
}

void Node::buildTree(int level) {
	auto bvh_sp = this->bvh.lock();
	switch (level % 3) {
	case 0:
		std::partial_sort(
			bvh_sp->triangles.begin() + this->from,
			bvh_sp->triangles.begin() + this->to,
			bvh_sp->triangles.end(),
			[level](Triangle* a, Triangle* b) -> bool { return a->getCenter().x > b->getCenter().x; }
		);
		break;
	case 1:
		std::partial_sort(
			bvh_sp->triangles.begin() + this->from,
			bvh_sp->triangles.begin() + this->to,
			bvh_sp->triangles.end(),
			[level](Triangle* a, Triangle* b) -> bool { return a->getCenter().y > b->getCenter().y; }
		);
		break;
	case 2:
		std::partial_sort(
			bvh_sp->triangles.begin() + this->from,
			bvh_sp->triangles.begin() + this->to,
			bvh_sp->triangles.end(),
			[level](Triangle* a, Triangle* b) -> bool { return a->getCenter().z > b->getCenter().z; }
		);
		break;
	}

	if (to - from < 5) {
		this->isLeaf_ = true;
		return;
	}

	std::vector<Division> divs;
	int avgLength = (this->to - this->from) / 5;
	int pivotIndex = this->from + avgLength;

	for (int i = 0; i < 4; i++) {
		divs.push_back(Division{
			std::make_shared<Node>(this->from, pivotIndex, this->bvh.lock()),
			std::make_shared<Node>(pivotIndex + 1, this->to, this->bvh.lock())
		});
		pivotIndex += avgLength;
	}

	Division bestDiv = divs[0];
	float bestDivVal = FLT_MAX;
	for (auto div : divs) {
		float divVal = calculateSAH(div);
		if (divVal < bestDivVal) {
			bestDivVal = divVal;
			bestDiv = div;
		}
	}

	this->children[0] = bestDiv.left;
	this->children[1] = bestDiv.right;

	this->children[0]->buildTree(level + 1);
	this->children[1]->buildTree(level + 1);
}

bool Node::isIntersected(RTCRay ray) {
	return this->getAABB()->isIntersected(ray);
}

std::shared_ptr<Triangle> Node::traverse(RTCRay & ray) {
	if (!this->isIntersected(ray)) {
		ray.tfar = FLT_MAX;
		return nullptr;
	}

	const Vector3 origin{ ray.org_x, ray.org_y, ray.org_z };
	auto bvh_sp = this->bvh.lock();
	
	if (this->isLeaf()) {
		std::shared_ptr<Triangle> bestTriangle = nullptr;
		float bestTriangleDist = FLT_MAX;

		for (int i = from; i <= this->to; i++) {
			auto hitPoint = bvh_sp->triangles[i]->intersect(ray);
			if (hitPoint != nullptr) {
				float dist = (origin - *hitPoint).L2Norm();
				if (dist < bestTriangleDist) {
					bestTriangleDist = dist;
					bestTriangle = std::make_shared<Triangle>(*bvh_sp->triangles[i]);
				}
			}
		}
		
		//if (bestTriangleDist != FLT_MAX) 
	//	std::cout << bestTriangleDist << std::endl;
		ray.tfar = bestTriangleDist;
		return bestTriangle;
	}

	RTCRay lRay(ray);
	RTCRay rRay(ray);
	auto lRes = this->children[0]->traverse(lRay);
	auto rRes = this->children[1]->traverse(rRay);

	if (lRay.tfar < rRay.tfar) {
		ray.tfar = lRay.tfar;
		return lRes;
	}
	ray.tfar = rRay.tfar;
	return rRes;
}

float Node::getSAHNodeVal(float parentSA) {
	const float t_traversal = 1.0f;
	const float t_intersect = 2.0f;

	auto bvh_sp = this->bvh.lock();
	float sum = 0.0f;
	for (int i = this->from; i <= this->to; i++) {
		sum += t_intersect * bvh_sp->triangles[i]->getArea();
	}
	return this->getHitProbability(parentSA) * sum;
}

float Node::getHitProbability(float parentSA) {
	return this->getAABB()->getSurfaceArea() / parentSA;
}






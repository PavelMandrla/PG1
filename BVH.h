#pragma once

#include <array>
#include <vector>
#include <algorithm>
#include <memory>
#include "triangle.h"

#include <iostream>

struct TriangleAndDistance {
	std::shared_ptr<Triangle> triangle;
	float distance;
};

struct AABB {
	Vector3 bounds[2];
	
	AABB(){}
	AABB(Vector3 v0, Vector3 v1) {
		bounds[0] = v0;
		bounds[1] = v1;
	}
	//float SurfaceArea() // for SAH
};

class BVH : public std::enable_shared_from_this<BVH> {
public:
	BVH(std::vector<std::shared_ptr<Triangle>> items);
	~BVH();

	void initialize();
	std::shared_ptr<Triangle> traverse(RTCRay & ray);

private:
	class Node {
	private:
		std::weak_ptr<BVH> bvh;
		std::shared_ptr<AABB> bounds;
	public:
		int span[2];
		std::shared_ptr<Node> children[2];

		Node(int from, int to, std::shared_ptr<BVH> bvh);

		std::shared_ptr<AABB> getAABB();
		TriangleAndDistance intersect(RTCRay &ray);

		bool isIntersected(RTCRay &ray);
		bool isLeaf();
	};
	std::shared_ptr<Node> BuildTree(int from, int to, int depth);
	
	std::shared_ptr<Node> root;
	std::vector<std::shared_ptr<Triangle>> items;
};


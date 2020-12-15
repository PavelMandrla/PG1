#pragma once

#include <vector>
#include <memory>
#include "triangle.h"
#include <embree3/rtcore_ray.h>

class AABB;
class Node;

struct Division {
	std::shared_ptr<Node> left;
	std::shared_ptr<Node> right;
};

class BVH : public std::enable_shared_from_this<BVH> {
public:
	BVH(std::vector<Triangle *> triangles);
	void initialize();
	std::shared_ptr<Triangle> traverse(RTCRay & ray);

	std::vector<Triangle *> triangles;
private:
	
	std::shared_ptr<Node> root;

	friend class Node;
	friend class AABB;
};

class Node {
public:
	Node(int from, int to, std::shared_ptr<BVH> bvh);
	bool isLeaf();
	std::shared_ptr<AABB> getAABB();
	
	void buildTree(int level);

	std::shared_ptr<Triangle> traverse(RTCRay & ray);
private:
	std::weak_ptr<BVH> bvh;
	std::shared_ptr<AABB> aabb;
	std::shared_ptr<Node> children[2];

	int from;
	int to;
	bool isLeaf_;

	float getSAHNodeVal(float parentSA);
	float getHitProbability(float parentSA);
	float calculateSAH(Division div);

	friend class AABB;
};

class AABB {
public:
	AABB(Node* node);
	bool isIntersected(RTCRay ray);
	float getSurfaceArea();
private:
	Vector3 bounds[2];
};

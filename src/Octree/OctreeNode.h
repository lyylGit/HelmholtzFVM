
//https://github.com/bertaye/Octree

#pragma once
#include <iostream>
#include <vector>
#include <memory>

using namespace std;
const int CHILD_COUNT = 8;

enum class OctreeNodeIndex
{
	/*
	*	______________
	*	|	   |	  |
		| LUB  | RUB  |
		|______|______|		LUB: Left Up Back
		|	   |	  |
		| LDB  | RDB  |
		|______|______|

		000 -> LEFT, DOWN, BACK
		111 -> RIGHT, UP, FRONT
	*/

	LeftDownBack = 0,
	LeftUpBack = 2,
	RightDownBack = 4,
	RightUpBack = 6,

	LeftDownFront = 1,
	LeftUpFront = 3,
	RightDownFront = 5,
	RightUpFront = 7,
};

template <class T>
class OctreeNode
{
private:
	OctreeNode *parent;
	float size;
	vector3D<T> position; //x,y,z poisiton of node 

    bool isLeaf;

	bool subNodesHasObject = false;
	unique_ptr <OctreeNode[]> subNodes;//[CHILD_COUNT]; // subnodes
	std::vector<long> objects;		   // objects of this node. this can be points, vertices, etc.
	
	
public:
	OctreeNode() : parent(nullptr),  size(-1), position(0, 0, 0),isLeaf(true) 
	{
	}

private:
	void SubdivideNode( float threshold)
	{
		isLeaf = false;
		float newPos[3] = {0, 0, 0};
		newPos[0] = position.x;
		newPos[1] = position.y;
		newPos[2] = position.z;

		subNodes=make_unique <OctreeNode[]>(CHILD_COUNT);

		for (int i = 0; i < 8; i++)
		{
			newPos[0] = position.x;
			newPos[1] = position.y;
			newPos[2] = position.z;

			// right
			if (i == 4 || i == 6 || i == 5 || i == 7)
				newPos[0] += size * 0.25f; 
			else
				newPos[0] -= size * 0.25f;
            // up
			if (i == 2 || i == 6 || i == 3 || i == 7)
				newPos[1] += size * 0.25f;
			else
				newPos[1] -= size * 0.25f;
            // front
			if (i == 1 || i == 3 || i == 5 || i == 7)
				newPos[2] += size * 0.25f;
			else
				newPos[2] -= size * 0.25f;

			subNodes[i].size = size * 0.5f;
			subNodes[i].position.x = newPos[0];
			subNodes[i].position.y = newPos[1];
			subNodes[i].position.z = newPos[2];

			subNodes[i].parent = this;
		}
	}
	
	int GetIndexByPosition(vector3D<T>&lookUpPosition) const
	{
		if (isLeaf)
		{
			// COMMANTABLE_OUTPUT("Node leaf cant give index by position!\n");
			return -1;
		}
		int index = 0;

		if (lookUpPosition.x> position.x)
			index += 4;
		if (lookUpPosition.y > position.y)
			index += 2;
		if (lookUpPosition.z > position.z)
			index += 1;
		return index;
	}
	void RecursivelyMarkParent()
	{
		if (parent != nullptr)
		{
			this->subNodesHasObject = true;
			parent->RecursivelyMarkParent();
		}
	}

	OctreeNode<T> *GetClosestNonEmptyNodeByPosition(vector3D<T> lookUpPosition) 
	{
		/*
		If we have a concave shape, then we cant just simply look at position by indexing. Because this may mislead
		Instead we need to find the shortest distanced and with object subnode of each node.
		We start with root, then we will compare all results coming from 8 child of root.
		Thus we will obtain the smallest node who is closest to the lookUpPosition.
	*/
		float tempDist[CHILD_COUNT];
		OctreeNode *tempNode[CHILD_COUNT];

		for (int i = 0; i < CHILD_COUNT; i++)
		{
			tempDist[i] = INFINITY;
			tempNode[i] = RecursiveNonEmptyFinder(lookUpPosition, &subNodes[i], tempDist[i]);
		}
		int shortestIdx = -1;
		float temp = INFINITY;
		for (int i = 0; i < CHILD_COUNT; i++)
		{
			if (temp > tempDist[i])
			{
				temp = tempDist[i];
				shortestIdx = i;
			}
		}

		return tempNode[shortestIdx];
	}
	OctreeNode<T> *RecursiveNonEmptyFinder(vector3D<T> &lookUpPosition, OctreeNode *node, float &lastDist)
	{
		if (node->isLeaf)
			return node;
		float currentDist = 0.0f, shortestDist = 1e9;
		OctreeNode *tempNode = node;
		for (int i = 0; i < CHILD_COUNT; i++)
		{
			if (node->subNodes[i].subNodesHasObject)
			{
				currentDist = sqrDistance(node->subNodes[i].position, lookUpPosition);
				if (shortestDist > currentDist)
				{
					shortestDist = currentDist;
					lastDist = shortestDist;
					tempNode = &node->subNodes[i];
				}
			}
		}
		if (tempNode != node)
			tempNode = RecursiveNonEmptyFinder(lookUpPosition, tempNode, lastDist);

		return tempNode;
	}
	float sqrDistance(const vector3D<T>& point1, const vector3D<T> &point2) const
	{
		float temp = (pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2) + pow(point1.z - point2.z, 2));
		return temp;
	}
	template<class U> friend class Octree;
};
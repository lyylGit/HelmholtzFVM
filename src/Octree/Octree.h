
//https://github.com/bertaye/Octree

#pragma once

#include "OctreeNode.h"
#include <typeinfo>
#include <math.h>
#include <memory>

using namespace std;

template <class T>
class Octree
{
private:
	float size; // size of biggest box
	float threshold;

	// 点集的边界坐标
	float leftBoundary, rigthBoundary,
		upBoundary, downBoundary,
		frontBoundary, backBoundary;

public:
	shared_ptr<vector3D<T>[]> points;
	long nCounts;

	vector3D<T> centerPoint;
	vector3D<T> geometricCenterPoint;

	OctreeNode<T> root;
	int totalSubdiv = 0;

public:
	Octree(shared_ptr<vector3D<T>[]> ptlist, long n) : size(-1), threshold(INFINITY), points(ptlist), nCounts(n)
	{
		calcuGeometric();
	}
	void calcuGeometric()
	{
		if (nCounts <= 1)
			return;

		int leftBoundaryPoint = 0,
			rigthBoundaryPoint = 0,
			upBoundaryPoint = 0,
			downBoundaryPoint = 0,
			frontBoundaryPoint = 0,
			backBoundaryPoint = 0;

		centerPoint = points[0];
		for (int i = 1; i < nCounts; i++)
		{
			auto p = points[i];
			centerPoint += p;

			if (p.x < points[leftBoundaryPoint].x)
				leftBoundaryPoint = i;
			if (p.x > points[rigthBoundaryPoint].x)
				rigthBoundaryPoint = i;
			if (p.y < points[downBoundaryPoint].y)
				downBoundaryPoint = i;
			if (p.y > points[upBoundaryPoint].y)
				upBoundaryPoint = i;
			if (p.z < points[backBoundaryPoint].z)
				backBoundaryPoint = i;
			if (p.z > points[frontBoundaryPoint].z)
				frontBoundaryPoint = i;
		}
		centerPoint /= (nCounts + 0.0);

		geometricCenterPoint.x = (points[leftBoundaryPoint].x + points[rigthBoundaryPoint].x) / 2.0f;
		geometricCenterPoint.y = (points[upBoundaryPoint].y + points[downBoundaryPoint].y) / 2.0f;
		geometricCenterPoint.z = (points[backBoundaryPoint].z + points[frontBoundaryPoint].z) / 2.0f;

		leftBoundary = points[leftBoundaryPoint].x;
		rigthBoundary = points[rigthBoundaryPoint].x;
		downBoundary = points[downBoundaryPoint].y;
		upBoundary = points[upBoundaryPoint].y;
		backBoundary = points[backBoundaryPoint].z;
		frontBoundary = points[frontBoundaryPoint].z;
	}
	float GetSize()
	{
		float tempSize = 0;

		tempSize = rigthBoundary - leftBoundary;

		if (upBoundary - downBoundary > tempSize)
			tempSize = upBoundary - downBoundary;
		if (frontBoundary - backBoundary > tempSize)
			tempSize = frontBoundary - backBoundary;

		return tempSize;
	}
	void BuiltOctree(float _threshold)
	{
		size = GetSize();
		threshold = _threshold;
		root.position = geometricCenterPoint;
		root.size = size;

		OctreeNode<T> *currentNode = &root;

		for (int i = 0; i < nCounts; i++)
		{
			if (sqrDistance(currentNode->position, points[i]) > threshold)
			{
				RecursiveSubdivision(currentNode, i);
				currentNode = &root; // we need to get back to root for new point.
			}
		}
	}
	void RecursiveSubdivision(OctreeNode<T> *_currentNode, long ptIndex)
	{
		if (_currentNode->size > threshold)
		{
			if (_currentNode->isLeaf)
			{
				_currentNode->SubdivideNode(threshold);
				totalSubdiv++;
			}
			_currentNode = &_currentNode->subNodes[_currentNode->GetIndexByPosition(points[ptIndex])];
			RecursiveSubdivision(_currentNode, ptIndex);
		}
		else
		{
			_currentNode->RecursivelyMarkParent();
			_currentNode->isLeaf = true;
			_currentNode->objects.push_back(ptIndex);
		}
		return;
	}
	float sqrDistance(vector3D<T> &point1, vector3D<T> &point2)
	{
		float temp = (pow(point1.x - point2.x, 2) + pow(point1.y - point2.y, 2) + pow(point1.z - point2.z, 2));
		return temp;
	}
	long GetClosestObject(T x,T y,T z)
	{
		vector3D<T> point(x,y,z);
		return GetClosestObject(point);
	}
	long GetClosestObject(vector3D<T> &point)
	{
		auto temp = root.GetClosestNonEmptyNodeByPosition(point);

		//寻找最近的点
		long closest = -1;
		float dist = 1e30;
		for (int i = 0; i < temp->objects.size(); i++)
		{
			auto pos = temp->objects[i];
			auto d2 = (point - points[pos]).norm2();
			if (d2 < dist)
			{
				closest = pos;
				dist = d2;
			}
		}

		return closest;
	}
};
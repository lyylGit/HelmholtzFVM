#pragma once

#include "Octree.h"

template <class TCase>
class MeshOctree
{
public:
    shared_ptr<Octree<TCase>> octree;
    // 构造器
public:
    MeshOctree(Mesh3<TCase> *_pMesh)
    {
        // 构造
        octree = make_shared<Octree<TCase>>(_pMesh->cellCenter, _pMesh->nCells);
        octree->BuiltOctree(0.010f);
    }

    long GetClosetCell(TCase x, TCase y, TCase z)
    {
        auto p = octree->GetClosestObject(x, y, z);
        return p;
    }
     long GetClosetCell(vector3D<TCase>& pos)
    {
        auto p = octree->GetClosestObject(pos.x, pos.y, pos.z);
        return p;
    }
};

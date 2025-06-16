#pragma once

#include "vector3d.h"
#include "CalcuGeo.h"
// #include "map"

using namespace std;

template <class T>
class Mesh3;

// extern double m_NearZero ;
template <class T>
class Cell
{
public:
    Cell(/* args */) {};
    ~Cell() {
    };
    // 面元集合
    vector<long> faces;
    // 面元法向的方向:1,-1
    vector<int> normalDir;
    // 相邻单元几何
    vector<long> neighbours;

    // 网格中心
    vector3D<T> _Center;
    // 体积
    T volume;
    // 网格中心到各个面的距离
    map<long, T> _distanceCenterToFaces;

public:
    // 添加单元信息
    bool AddFace(long faceID, int normaldir)
    {
        if (find(faces.begin(), faces.end(), faceID) != faces.end())
        {
            throw "该面元已经存在,不能重复添加";
            return false;
        }

        faces.push_back(faceID);
        normalDir.push_back(normaldir);
        return true;
    }
    
    // 计算体积和各个面的距离
    void CalcuVolume(Mesh3<T> *pMesh, int iCell)
    {
        // 先把 faces中的点找齐
        vector<vector<long>> _faces;
        int nfaces = faces.size();
        for (int j = 0; j < nfaces; j++)
        {
            vector<long> face;
            long i = faces[j];
            // 该面点数
            auto npts = pMesh->faceStructrue[i];
            // 该面的结构
            for (long k = 0; k < npts; k++)
                face.push_back(pMesh->faceIndex[i][k]);
            // 需要调整为边界法向向外
            if (pMesh->nextCell[i] == -1)
                reverse(face.begin(), face.end());

            _faces.push_back(face);
        }
        // 取个顶点
        long iTopPoint = _faces[0][0];
        // 顶点坐标
        // VectorMY<T> xyzTopPoint(pMesh->nodes, iTopPoint);
        vector3D<T> xyzTopPoint(pMesh->nodes, iTopPoint);
        volume = 0;
        // 同时计算中心
        _Center.Zero();
        _Center += xyzTopPoint;
        int n = 1;
        for (int i = 1; i < nfaces; i++)
        {
            long curFaceId = faces[i];
            auto f = _faces[i];
            if (find(f.begin(), f.end(), iTopPoint) == f.end())
            {
                auto points = pMesh->FacePoints(f);
                volume += PyramidVolume(points, xyzTopPoint);
                int nn = points.size();
                for (int k = 0; k < nn; k++)
                {
                    _Center += points[k];
                    n++;
                }
                // Delete(points);
            }
            // 计算相邻单元
            if (pMesh->prevCell[curFaceId] * pMesh->nextCell[curFaceId] >= 0)
            {
                if (pMesh->prevCell[curFaceId] == iCell)
                    neighbours.push_back(pMesh->nextCell[curFaceId]);
                else
                    neighbours.push_back(pMesh->prevCell[curFaceId]);
            }
        }

        _Center = _Center / n;

        // 计算中点到各个面的距离
        for (size_t i = 0; i < faces.size(); i++)
        {
            long curFaceId = faces[i];
            auto f = _faces[i];
            auto points = pMesh->FacePoints(f);
            _distanceCenterToFaces[curFaceId] = PyramidHeight(points, _Center);
            // Delete(points);
        }
    }
};
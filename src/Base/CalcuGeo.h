#pragma once
#include <vector>
#include "vector3d.h"
#include "Types.h"
using namespace std;

template <class TCase, class TField>
class FieldBase;

template <class TCase>
class Mesh3;

class EquationConfigue;

extern double m_NearZero;
/*
template <class T>
vector<vector3D<T>> Copy(vector<vector3D<T>> &array);

template <class T>
vector3D<T> CalcuNormal(vector<vector3D<T>> &saPoints);

template <class T>
vector<vector3D<T>> LocalCoordinates(vector<vector3D<T>> &saPoints);

template <class T>
T PolygonArea(vector<vector3D<T>> saPoints);

template <class T>
T PyramidHeight(vector<vector3D<T>> &points, vector3D<T> &point);

template <class T>
T PyramidVolume(vector<vector3D<T>> &points, vector3D<T> &point);
*/

template <class T>
vector<vector3D<T>> Copy(vector<vector3D<T>> &array)
{
    vector<vector3D<T>> ret;
    int n = array.size();
    for (int i = 0; i < n; i++)
        ret.push_back(vector3D<T>(array[i].x, array[i].y, array[i].z));

    return ret;
}
template <class T>
vector3D<T> CalcuCenter(vector<vector3D<T>> &saPoints)
{
    vector3D<T> center(0, 0, 0);
    int n = saPoints.size();
    for (int i = 0; i < n; i++)
        center += saPoints[i];

    center = center / n;

    return center;
}
template <class T>
vector3D<T> CalcuNormal(vector<vector3D<T>> &saPoints)
{
    vector3D<T> zero(0, 0, 0);
    if (saPoints.size() < 3)
    {
        return zero;
    }
    if (saPoints.size() == 3) // 三个点
    {
        vector3D<T> normal, v1, v2;
        v1 = (saPoints[1]) - (saPoints[0]);
        v2 = (saPoints[2]) - (saPoints[1]);
        normal = v1.cross(v2);
        normal.unit();
        if (normal.norm() < m_NearZero)
            return zero;
        return normal;
    }
    else // 取最佳平面
    {
        vector3D<T> normal(0, 0, 0);
        auto p = saPoints[saPoints.size() - 1];
        int n = saPoints.size();
        for (int i = 0; i < n; ++i)
        {
            auto c = saPoints[i];
            normal.x += (p.z + c.z) * (p.y - c.y);
            normal.y += (p.x + c.x) * (p.z - c.z);
            normal.z += (p.y + c.y) * (p.x - c.x);
            p = c;
        }
        normal.unit();

        if (normal.norm() < m_NearZero)
            return zero;
        return normal;
    }
    return zero;
}

template <class T>
vector<vector3D<T>> LocalCoordinates(vector<vector3D<T>> &saPoints)
{
    auto points = Copy(saPoints);

    // Select p.points[0] as the displacement
    auto Rx = points[0].x;
    auto Ry = points[0].y;
    auto Rz = points[0].z;

    // Subtract R from all the points of polygon p
    int n = points.size();
    for (int k = 0; k < n; k++)
    {
        points[k].x -= Rx;
        points[k].y -= Ry;
        points[k].z -= Rz;
    }

    // Select P0P1 as the x-direction
    vector3D<T> iprime;
    iprime = points[1] - points[0];

    // Find a unit vector in the xprime direction
    iprime.unit();

    // Find the vector P1P2
    vector3D<T> p1p2;
    p1p2 = points[2] - points[1];

    // Find a vector kprime in the zprime direction
    vector3D<T> kprime;
    kprime = iprime.cross(p1p2);

    // Make kprime a unitVector
    kprime.unit();

    // Find the vector jprime in the yprime direction
    vector3D<T> jprime;
    jprime = kprime.cross(iprime);

    // For each point, calculate the projections on xprime, yprime, zprime
    // (All zprime values should be zero)
    for (int k = 0; k < n; k++)
    {
        // VectorMY<T> pprime(0, 0, 0);
        vector3D<T> pv(points[k].x, points[k].y, points[k].z);
        points[k].x = iprime * pv;
        points[k].y = jprime * pv;
        points[k].z = kprime * pv;
    }

    // Return a polygon in its own local x'y'z' coordinates
    return points;
}

template <class T>
T PolygonArea(vector<vector3D<T>> saPoints)
{
    // Get the polygon in its own
    // x-y coordinates (all z's should be 0)

    auto points = LocalCoordinates(saPoints);

    // Apply the surveyor's formula
    int len = points.size();

    T result = 0.0;
    T dx = 0.0;
    T dy = 0.0;
    for (int k = 0; k < (len - 1); k++)
    {
        dx = points[k + 1].x - points[k].x;
        dy = points[k + 1].y - points[k].y;
        result += points[k].x * dy -
                  points[k].y * dx;
    }
    dx = points[0].x - points[len - 1].x;
    dy = points[0].y - points[len - 1].y;
    result += points[len - 1].x * dy -
              points[len - 1].y * dx;

    // Delete(points);
    return result / 2.0;
}
template <class T>
T PyramidHeight(vector<vector3D<T>> &points, vector3D<T> &point)
{
    // Construct a vector from the Pyramid base to the top point
    vector3D<T> vt;
    vt = point - points[0];

    // Calculate the perpendicular
    // distance from the base to the top point.
    // The distance d is the projection of vt in the normal direction.
    // Because a right-hand coordinate system is assumed, the value of d
    // may be negative, so the absolute value is returned.
    vector3D<T> norm;
    norm = CalcuNormal(points);
    auto d = norm * vt;
    auto result = 0.0;
    if (d < 0.0)
        result = fabs(d);
    else
        result = d;
    return result;
}

template <class T>
T PyramidVolume(vector<vector3D<T>> &points, vector3D<T> &point)
{
    // Calculate the perpendicular distance
    // from the base to the top point
    auto d = PyramidHeight(points, point);

    // Calculate the area of the base
    auto baseArea = PolygonArea(points);

    // Calculate the volume of the polygon's pyramid
    auto volume = d * baseArea / 3.0;
    return volume;
}

// 计算交界面值，由于多处用到，在此定义
template <class TCase, class TField>
bool CalcuFaceValues(Mesh3<TCase> *_Mesh, EquationConfigue *_EQConfig, const FieldBase<TCase, TField> *_Field)
{
    size_t n = _Mesh->inFaces.size();
    // 中心格式
    if (_EQConfig == NULL || _EQConfig->dmsFaceValue == DiscreteMethodSchemeType::Central)
    {
        for (size_t i = 0; i < n; i++)
        {
            long faceId = _Mesh->inFaces[i];
            auto fC = _Mesh->ownerCell[faceId];
            auto fF = _Mesh->neighborCell[faceId];
            auto gc = _Mesh->gc[faceId];
            (*_Field->fFaceValues)[faceId] = gc * (*_Field->fValues)[fC] + (1 - gc) * (*_Field->fValues)[fF];
        }
        return true;
    }
    else
    {
        throw "计算交界面上的值，其他格式尚未支持";
        return false;
    }
}
//矢量
template <class TCase, class TField>
bool CalcuFaceValuesVec(Mesh3<TCase> *_Mesh, EquationConfigue *_EQConfig, const FieldBase<TCase, vector3D<TField>> *_Field)
{
    size_t n = _Mesh->inFaces.size();
    // 中心格式
    if (_EQConfig == NULL || _EQConfig->dmsFaceValue == DiscreteMethodSchemeType::Central)
    {
        for (size_t i = 0; i < n; i++)
        {
            long faceId = _Mesh->inFaces[i];
            auto fC = _Mesh->ownerCell[faceId];
            auto fF = _Mesh->neighborCell[faceId];
            auto gc = _Mesh->gc[faceId];
            (*_Field->fFaceValues)[faceId] = gc * (*_Field->fValues)[fC] + (1 - gc) * (*_Field->fValues)[fF];
        }
        return true;
    }
    else
    {
        throw "计算交界面上的值，其他格式尚未支持";
        return false;
    }
}
// 从控制体梯度GradV去计算交界面梯度GradF
// 处理边界上的梯度
template <class TCase, class TField, class TOther>
void _InterpolateGradV2F_BC(const FieldBase<TCase, TOther> *_Field, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF)
{
    auto _Case = _Field->_Case;
    auto mesh = _Case->_Mesh;
    for (auto _bc : _Case->PhysicsBCs)
    {
        auto bc = dynamic_pointer_cast<PhysicsBCBase>(_bc);
        long from = get<0>(bc->_BCFaces), to = get<1>(bc->_BCFaces);
        if (bc->_BCType != PhysicsBCType::Symmetry)
        {
            for (long faceId = from; faceId <= to; faceId++)
            {
                auto fC = mesh->ownerCell[faceId];
                auto Ecf = mesh->ECF[faceId];
                auto dcf = mesh->dcf[faceId];

                // gradF[faceId] = -GradV[fC].dot(Ecf) * Ecf + ((*_Field->fFaceValues)[faceId] - (*_Field->fValues)[fC]) / dcf * Ecf;
                auto a = -GradV[fC].dot(Ecf) + ((*_Field->fFaceValues)[faceId] - (*_Field->fValues)[fC]) / dcf;
                gradF[faceId] = vector3D<TField>(a * Ecf.x, a * Ecf.y, a * Ecf.z);
            }
        }
    }
}
template <class TCase, class TField, class TOther>
void _InterpolateGradV2F(const FieldBase<TCase, TOther> *_Field, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF, FaceValueSchemeType _AvgType)
{
    auto _Case = _Field->_Case;
    auto mesh = _Case->_Mesh;
    if (_AvgType == FaceValueSchemeType::Central)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];
        }
    }
    else if (_AvgType == FaceValueSchemeType::Corrected)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];

            auto Ecf = mesh->ECF[faceId];
            auto dcf = mesh->dcf[faceId];
            auto local = ((*_Field->fValues)[fF] - (*_Field->fValues)[fC]) / dcf * Ecf;
            // auto avg = gradF[faceId].dot(Ecf) * Ecf;
            auto _avg = gradF[faceId].dot(Ecf);
            vector3D<TField> avg(Ecf.x * _avg, Ecf.y * _avg, Ecf.z * _avg);
            gradF[faceId] += -avg + local;
        }
    }

    _InterpolateGradV2F_BC(_Field, GradV, gradF);
}
template <class TCase, class TField, class TOther>
void _InterpolateGradV2F_BC(const TOther v, CaseBase<TCase> *_Case, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF)
{
    // auto _Case = _Field->_Case;
    auto mesh = _Case->_Mesh;
    for (auto _bc : _Case->PhysicsBCs)
    {
        auto bc = dynamic_pointer_cast<PhysicsBCBase>(_bc);
        long from = get<0>(bc->_BCFaces), to = get<1>(bc->_BCFaces);
        if (bc->_BCType != PhysicsBCType::Symmetry)
        {
            for (long faceId = from; faceId <= to; faceId++)
            {
                auto fC = mesh->ownerCell[faceId];
                if (fC < 0)
                    fC = mesh->neighborCell[faceId];
                auto Ecf = mesh->ECF[faceId];
                auto dcf = mesh->dcf[faceId];

                // gradF[faceId] = -GradV[fC].dot(Ecf) * Ecf + ((*_Field->fFaceValues)[faceId] - (*_Field->fValues)[fC]) / dcf * Ecf;
                auto a = -GradV[fC].dot(Ecf); //+ ((*_Field->fFaceValues)[faceId] - (*_Field->fValues)[fC]) / dcf;
                gradF[faceId] = vector3D<TField>(a * Ecf.x, a * Ecf.y, a * Ecf.z);
            }
        }
    }
}
template <class TCase, class TField, class TOther>
void _InterpolateGradV2F(const TOther v, CaseBase<TCase> *_Case, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF, FaceValueSchemeType _AvgType)
{
    auto mesh = _Case->_Mesh;
    if (_AvgType == FaceValueSchemeType::Central)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];
        }
    }
    else if (_AvgType == FaceValueSchemeType::Corrected)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];

            auto Ecf = mesh->ECF[faceId];
            auto dcf = mesh->dcf[faceId];
            // TField local =0;// ((*_Field->fValues)[fF] - (*_Field->fValues)[fC]) / dcf * Ecf;
            // auto avg = gradF[faceId].dot(Ecf) * Ecf;
            auto _avg = gradF[faceId].dot(Ecf);
            vector3D<TField> avg(Ecf.x * _avg, Ecf.y * _avg, Ecf.z * _avg);
            gradF[faceId] += -avg; // + local;
        }
    }

    _InterpolateGradV2F_BC(v, _Case, GradV, gradF);
}

template <class TCase, class TField, class TOther>
void _InterpolateGradV2F_BC(const TOther *_Field, CaseBase<TCase> *_Case, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF)
{
    // auto _Case = _Field->_Case;
    auto mesh = _Case->pMesh;
    for (auto _bc : _Case->PhysicsBCs)
    {
        auto bc = dynamic_pointer_cast<PhysicsBCBase>(_bc);
        long from = get<0>(bc->_BCFaces), to = get<1>(bc->_BCFaces);
        if (bc->_BCType != PhysicsBCType::Symmetry)
        {
            for (long faceId = from; faceId <= to; faceId++)
            {
                auto fC = mesh->ownerCell[faceId];
                auto Ecf = mesh->ECF[faceId];
                auto dcf = mesh->dcf[faceId];

                // gradF[faceId] = -GradV[fC].dot(Ecf) * Ecf + ((*_Field->fFaceValues)[faceId] - (*_Field->fValues)[fC]) / dcf * Ecf;
                auto a = -GradV[fC].dot(Ecf) + (_Field[faceId] - _Field[fC]) / dcf;
                gradF[faceId] = vector3D<TField>(a * Ecf.x, a * Ecf.y, a * Ecf.z);
            }
        }
    }
}
template <class TCase, class TField, class TOther>
void _InterpolateGradV2F(const TOther *_Field, CaseBase<TCase> *_Case, shared_ptr<vector3D<TField>[]> GradV, shared_ptr<vector3D<TField>[]> gradF, FaceValueSchemeType _AvgType)
{
    // auto _Case = _Field->_Case;
    auto mesh = _Case->pMesh;
    if (_AvgType == FaceValueSchemeType::Central)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];
        }
    }
    else if (_AvgType == FaceValueSchemeType::Corrected)
    {
        for (long faceId : mesh->inFaces)
        {
            auto fC = mesh->ownerCell[faceId];
            auto fF = mesh->neighborCell[faceId];
            auto gc = mesh->gc[faceId];
            gradF[faceId] = gc * GradV[fC] + (1 - gc) * GradV[fF];

            auto Ecf = mesh->ECF[faceId];
            auto dcf = mesh->dcf[faceId];
            auto local = (_Field[fF] - _Field[fC]) / dcf * Ecf;
            // auto avg = gradF[faceId].dot(Ecf) * Ecf;
            auto _avg = gradF[faceId].dot(Ecf);
            vector3D<TField> avg(Ecf.x * _avg, Ecf.y * _avg, Ecf.z * _avg);
            gradF[faceId] += -avg + local;
        }
    }

    _InterpolateGradV2F_BC(_Field, GradV, gradF);
}
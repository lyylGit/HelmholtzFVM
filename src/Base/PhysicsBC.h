#pragma once
#include <map>
#include "BC.h"

using namespace std;

class PhysicsBCBase
{
public:
    // 变量名称
    string Name;
    // 对应mesh->faceBD中的序号
    long ID;
    // 边界面集合:from,to,bcType,faceType;
    tuple<long, long, int, int> _BCFaces;
    // 边界类型
    PhysicsBCType _BCType;
    // 编号
};

template <class T>
class PhysicsBC : public PhysicsBCBase
{
public:
    // 边界静压
    double _pb;
    // 边界总压
    double _p0;
    // 边界流量
    double _mb;
    // 边界体积作用力
    vector3D<T> _Fb;
    // 边界速度
    vector3D<T> _Vb;
    // 边界速度方向
    vector3D<T> _Vbdir;

    // 其他标量边界值
    map<string, T> _htPhi;
};
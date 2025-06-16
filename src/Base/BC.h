#pragma once
#include <map>
#include "IName.h"
#include "MemSegmentBC.h"

using namespace std;

template <class TCase,class TField>
class Impedance;
template <class TCase>
class CaseBase;
template <class TCase, class TField>
class FieldBase;

// 对于不同模板参数的对象，采用名称作为唯一标志
class BCBase
{
public:
    // 变量名称
    string Name;
    // 对应mesh->faceBD中的序号
    long ID;
    // 索引信息
    long length, start, end;

    BoundaryConditionType _BCType;

    BCBase(string name, long id)
        : Name(name), ID(id)
    {
    }
    void InitInfo(tuple<long, long, int, int> &BCFaces)
    {
        start = get<0>(BCFaces), end = get<1>(BCFaces);
        length = end - start + 1;
    }
};

template <class TCase, class TField>
class BC : public BCBase
{
public:
    // 对应的场
    FieldBase<TCase, TField> *_Field;

    // 固定值/边界面法向梯度值:第1、2类边界条件
    // TField _SingleValue;
    shared_ptr<MemSegmentBC<TCase, TField>> _SingleValue12,
        // 混合边界第3类边界条件:A+Bϕ+Γ∇ϕ=0
        //_FieldA=0-->A=0;
        //_FieldB=0 _FieldGamma=0-->B=1 Γ=1;
        _FieldA3, _FieldB3, _FieldGamma3;
    // 定义边界条件上的声阻抗
    shared_ptr<Impedance<TCase,TField>> _Imp=nullptr;

public:
    BC(string name, long id)
        : BCBase(name, id)
    {
    }
    ~BC() {}

    void Allocate(shared_ptr<MemSegmentBC<TCase, TField>> &mem)
    {
        mem = make_shared<MemSegmentBC<TCase, TField>>(length, start);
    }
    void SetValue12(TField f)
    {
        Allocate(_SingleValue12);
        _SingleValue12->SetValues(f);
    }
    void SetValueA3(TField f)
    {
        Allocate(_FieldA3);
        _FieldA3->SetValues(f);
    }
    void SetValueB3(TField f)
    {
        Allocate(_FieldB3);
        _FieldB3->SetValues(f);
    }
    void SetValueGamma3(TField f)
    {
        Allocate(_FieldGamma3);
        _FieldGamma3->SetValues(f);
    }
    void SetValueGamma3(shared_ptr<Impedance<TCase,TField>> imp)
    {
        _Imp = imp;
        Allocate(_FieldGamma3);
        // _FieldGamma3->SetValues(f);
    }
    void CalcuImpValues(double fHz)
    {
        _Imp->CalcuImpValues(fHz, _Imp.get(), _FieldGamma3.get());
    }
    void Init()
    {
        if (this->_BCType == BoundaryConditionType::Value || this->_BCType == BoundaryConditionType::Gradient)
        {
            for (long i = start; i <= end; i++)
                (*_Field->fFaceValues)[i] = (*_SingleValue12)[i];
        }
    }
};

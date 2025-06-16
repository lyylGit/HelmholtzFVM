#pragma once

#include <cstddef>

using namespace std;

template <class TCase, class TField>
class FieldBase;

template <class TCase>
class Mesh3;

class EquationConfigue;

template <class TCase, class TField>
class CalcuBase
{
public:
    // 方程设置
    EquationConfigue *_EQConfig;
    // 计算的场
    const FieldBase<TCase, TField> *_Field;

public:
    CalcuBase(const FieldBase<TCase, TField> *pField, EquationConfigue *eqC)
        : _Field(pField), _EQConfig(eqC)
    {
    }

    ~CalcuBase() {}

public:
    // 公共函数
    //  算子：计算交界面上的值
    bool CalcuFaceValues()
    {
        return ::CalcuFaceValues(_Field->_Case->pMesh, _EQConfig, _Field);
    }
   
    // 继承类实现
    virtual bool Calcu() { return true; }
};
#pragma once

#include <cstddef>

template <class TCase, class TField>
class EquationBase;

template <class TCase, class TField>
class FieldBase;

template <class TCase>
class Mesh3;

class EquationConfigue;

template <class TCase, class TField>
class DiscreteBase
{

public:
    EquationBase<TCase, TField> *_Equation;

protected:
    shared_ptr<Flux<TCase, TField>> _Flux;

public:
    DiscreteBase(EquationBase<TCase, TField> *EQ)
    {
        _Equation = EQ;

        // 找到计算矩阵
        _Flux = EQ->_Flux;
    }
    ~DiscreteBase() {}

    // 采用中心格式计算交界面的值
    template <class TOther>
    bool CalcuFaceValues(const FieldBase<TCase, TOther> *pField)
    {
        return ::CalcuFaceValues<TOther>(_Equation->_Case->_Mesh, &_Equation->_EquationConfigue, pField);
    }
    template <class TOther>
    bool CalcuFaceValuesVec(const FieldBase<TCase, vector3D<TOther>> *pField)
    {
        return ::CalcuFaceValuesVec<TOther>(_Equation->_Case->_Mesh, &_Equation->_EquationConfigue, pField);
    }
    template <class TOther>
    void InterpolateGradV2F(const FieldBase<TCase, TOther> *pField, shared_ptr<vector3D<TField>[]> gradV, shared_ptr<vector3D<TField>[]> gradF, FaceValueSchemeType fs)
    {
        if (pField != 0)
            ::_InterpolateGradV2F<TCase, TField, TOther>(pField, gradV, gradF, fs);
        else
            ::_InterpolateGradV2F<TCase, TField, TOther>(1, _Equation->_Case, gradV, gradF, fs);
    }

    // 由子类实现
    // template <class TOther>
    // virtual bool Discrete(const FieldBase<TCase, TOther> *tagetField, const FieldBase<TCase, TOther> *gamma)
    // {
    //     return true;
    // };
};

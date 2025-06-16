#pragma once
using namespace std;
#include "PhysicsBC.h"

template <class TCase, class TField>
class DiscreteBase;

template <class TCase, class TField>
class FieldBase;

template <class TCase, class TField>
class DiscreteGrad : public DiscreteBase<TCase, TField>
{
public:
    // 非正交修正的限制子
    TCase lamda = 1 / 3.0;
    // 非正交修正的前一时刻梯度
    // 体积cell的梯度，面梯度
    shared_ptr<vector3D<TField>[]> gradV, gradF;

private:
    const FieldBase<TCase, TField> *_phiField;

public:
    DiscreteGrad(EquationBase<TCase, TField> *EQ)
        : DiscreteBase<TCase, TField>(EQ)
    {
    }
    ~DiscreteGrad() {}
    // Γ∙∇ϕ  :其中 Γ为送入的矢量场,ϕ为待求解场phiField
    // 将其转化为 sigma_f[(Γ_f∙S_f*ϕ_f)]
    // sigma_f 对网格cell的所有面f进行求和
    //_f为该面上的值，其中S_f为面积矢量
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const vector3D<TOther> &gamma)
    {
        // _gamma = gamma;
        _phiField = phiField;

        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        // 内部面
        Inner(gamma);

        // 边界面上的处理,
        for (long _bcId : _Mesh->FaceBD)
        {
            BC<TCase, TField> *bc = (BC<TCase, TField> *)_Case->template GetBoundaryConditions<TField>(phiField, _bcId);
            if (bc == 0)
                throw("边界条件未定义");
            else
                Boundary(bc, gamma);
        }

        return true;
    }
    template <class TOther>
    bool Inner(const vector3D<TOther> &gamma_f)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;

        for (auto faceId : _Mesh->inFaces)
        {
            auto fF = _Mesh->neighborCell[faceId];
            auto fC = _Mesh->ownerCell[faceId];
            auto gc = _Mesh->gc[faceId];

            auto S = _Mesh->surfaceVector[faceId];

            auto cf = gamma_f.dot(S);
            auto ff = -cf * (1 - gc);
            cf = cf * gc;

            (*flux->FluxC_f)[faceId] += cf;
            (*flux->FluxF_f)[faceId] += ff;
        }

        return true;
    }
    template <class TOther>
    bool Boundary(BC<TCase, TField> *bc, const vector3D<TOther> &gamma)
    {
        switch (bc->_BCType)
        {
        case BoundaryConditionType::Value:
            return Value(bc, gamma);
        case BoundaryConditionType::Gradient:
            return Gradient(bc, gamma);
        case BoundaryConditionType::Mixed:
            return Mixed(bc, gamma);
        case BoundaryConditionType::MixedModal:
            return MixedModal(bc, gamma);
        default:
            throw("该边界条件未实现");
        }
    }

    // 边界给定值
    template <class TOther>
    bool Value(BC<TCase, TField> *bc, const vector3D<TOther> &gamma_f)
    {
        // 边界给定值，对此无影响
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;

        for (auto faceId = from; faceId <= to; faceId++)
        {
            auto fC = _Mesh->ownerCell[faceId];
            auto S = _Mesh->surfaceVector[faceId];

            auto f = (*_phiField->fFaceValues)[faceId];
            auto cf = gamma_f.dot(S) * f;
            (*flux->FluxT_f)[faceId] += -cf;
        }

        return true;
    }
    // 边界给定梯度值
    template <class TOther>
    bool Gradient(BC<TCase, TField> *bc, const vector3D<TOther> &gamma_f)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;

        for (auto faceId = from; faceId <= to; faceId++)
        {
            auto fC = _Mesh->ownerCell[faceId];
            auto dcf = _Mesh->dcf[faceId];

            auto S = _Mesh->surfaceVector[faceId];
            // 边界是梯度
            auto f = (*_phiField->fFaceValues)[faceId];

            auto cf = gamma_f.dot(S);

            auto tf = f * dcf * cf;

            (*flux->FluxC_f)[faceId] += cf;
            (*flux->FluxT_f)[faceId] += -tf;
        }

        return true;
    }
    // 边界为第三类边界条件
    template <class TOther>
    bool Mixed(BC<TCase, TField> *bc, const vector3D<TOther> &gamma_f)
    {
        // throw("该边界条件未实现");
        return true;
    }
    // 边界为第三类边界条件计算模态时，只处理常数项，该项并入刚度矩阵
    template <class TOther>
    bool MixedModal(BC<TCase, TField> *bc, const vector3D<TOther> &gamma_f)
    {
        // throw("该边界条件未实现");
        return true;
    }
};
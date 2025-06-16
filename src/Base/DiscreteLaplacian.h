#pragma once
using namespace std;
#include "PhysicsBC.h"

template <class TCase, class TField>
class DiscreteBase;

template <class TCase, class TField>
class FieldBase;

template <class TCase, class TField>
class DiscreteLaplacian : public DiscreteBase<TCase, TField>
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
    DiscreteLaplacian(EquationBase<TCase, TField> *EQ)
        : DiscreteBase<TCase, TField>(EQ)
    {
    }
    ~DiscreteLaplacian() {}
    // ∇∙(Γ∇ϕ)  :其中 Γ为送入的场,ϕ为待求解场phiField
    // 将其转化为 sigma_f[(Γ_f*S_f*grad(ϕ)_f)]
    // sigma_f 对网格cell的所有面f进行求和
    //_f为该面上的值
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const FieldBase<TCase, TOther> *gamma)
    {
        // _gamma = gamma;
        _phiField = phiField;

        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        if (!this->_Equation->_Signal && gamma != 0)
        {
            // 计算Γ的界面值
            if (!this->template CalcuFaceValues<TOther>(gamma))
                return false;

            // 下次不需要计算
            this->_Equation->_Signal = true;
        }
        // 计算前一时刻的控制体梯度
        CalcuVolumeGrad<TCase, TField> cvg(phiField, &this->_Equation->_EquationConfigue);
        if (!cvg.Calcu())
            return false;

        gradV = cvg._Grad;
        gradF = make_shared<vector3D<TField>[]>(_Mesh->nFaces);
        this->template InterpolateGradV2F<TOther>(gamma, gradV, gradF, FaceValueSchemeType::Corrected);

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
    bool Inner(const FieldBase<TCase, TOther> *gamma)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (gamma != 0)
            for (auto faceId : _Mesh->inFaces)
            {
                auto fF = _Mesh->neighborCell[faceId];
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto E = _Mesh->ECF[faceId];
                auto dcf = _Mesh->dcf[faceId];
                // 非正交修正
                auto Ef = s * s / S.dot(E) * E;
                auto Tf = S - Ef;
                auto gdiff = s / dcf;

                auto gamm_f = (*gamma->fFaceValues)[faceId];
                auto cf = -sign * gamm_f * gdiff;
                auto ff = -cf;
                auto vf = sign * gamm_f * gradF[faceId] * Tf;
                auto tf = cf * (*_phiField->fValues)[fC] + ff * (*_phiField->fValues)[fF] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxF_f)[faceId] += ff;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }
        else
            for (auto faceId : _Mesh->inFaces)
            {
                auto fF = _Mesh->neighborCell[faceId];
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId]; // _Mesh->faceNormal[faceId] * s;
                auto E = _Mesh->ECF[faceId];
                auto dcf = _Mesh->dcf[faceId];
                // 非正交修正
                auto Ef = s * s / S.dot(E) * E;
                auto Tf = S - Ef;
                auto gdiff = s / dcf;

                TOther gamm_f = 1; // (*gamma->fFaceValues)[faceId];
                auto cf = -sign * gamm_f * gdiff;
                auto ff = -cf; // sign * gamm_f * gdiff;
                auto vf = sign * gamm_f * gradF[faceId] * Tf;
                auto tf = cf * (*_phiField->fValues)[fC] + ff * (*_phiField->fValues)[fF] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxF_f)[faceId] += ff;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }

        return true;
    }
    template <class TOther>
    bool Boundary(BC<TCase, TField> *bc, const FieldBase<TCase, TOther> *gamma)
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
    bool Value(BC<TCase, TField> *bc, const FieldBase<TCase, TOther> *gamma)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;
        TField sign = 1;
        if (gamma != 0)
            for (auto faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto E = _Mesh->ECF[faceId];
                auto dcf = _Mesh->dcf[faceId];
                // 非正交修正
                auto Ef = s * s / S.dot(E) * E;
                auto Tf = S - Ef;
                auto gdiff = s / dcf;

                auto gamm_f = (*gamma->fFaceValues)[faceId];
                auto cf = -sign * gamm_f * gdiff;
                auto ff = -cf * (*_phiField->fFaceValues)[faceId];
                auto vf = ff + sign * gamm_f * gradF[faceId] * Tf;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                //    (*flux->FluxF_f)[faceId] += 0;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }
        else
            for (auto faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto E = _Mesh->ECF[faceId];
                auto dcf = _Mesh->dcf[faceId];
                // 非正交修正
                auto Ef = s * s / S.dot(E) * E;
                auto Tf = S - Ef;
                auto gdiff = s / dcf;

                TOther gamm_f = 1;
                auto cf = -sign * gamm_f * gdiff;
                auto ff = -cf * (*_phiField->fFaceValues)[faceId];
                auto vf = sign * gamm_f * gradF[faceId] * Tf + ff;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxF_f)[faceId] += 0;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }

        return true;
    }
    // 边界给定梯度值
    template <class TOther>
    bool Gradient(BC<TCase, TField> *bc, const FieldBase<TCase, TOther> *gamma)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;
        TField sign = 1;
        if (gamma != 0)
            for (auto faceId = from; faceId <= to; faceId++)
            {

                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                // 边界是梯度
                auto f = (*_phiField->fFaceValues)[faceId];
                gradF[faceId] = vector3D<TField>(f * S.x, f * S.y, f * S.z);

                auto gamm_f = (*gamma->fFaceValues)[faceId];
                auto cf = -sign * gamm_f * gradF[faceId] * S;
                auto vf = cf;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }
        else
            for (auto faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];

                // 边界是梯度
                auto f = (*_phiField->fFaceValues)[faceId];
                gradF[faceId] = vector3D<TField>(f * S.x, f * S.y, f * S.z);

                TOther gamm_f = 1;
                auto cf = -sign * gamm_f * gradF[faceId] * S;
                auto vf = cf;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }

        return true;
    }
    // 边界为第三类边界条件
    template <class TOther>
    bool Mixed(BC<TCase, TField> *bc, const FieldBase<TCase, TOther> *gamma)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;
        TField sign = 1;
        if (gamma != 0)
            for (auto faceId = from; faceId <= to; faceId++)
            {

                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto dcf = _Mesh->dcf[faceId];

                auto gamm_f = (*gamma->fFaceValues)[faceId];
                // 边界定义
                auto A = (*bc->_FieldA3)[faceId];
                auto B = (*bc->_FieldB3)[faceId];
                auto C = (*bc->_FieldGamma3)[faceId];

                auto K = C / dcf;
                auto c = -sign * gamm_f * K / (C * (B + K)) * s;
                auto cf = B * c;
                auto vf = A * c;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }
        else
            for (auto faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto dcf = _Mesh->dcf[faceId];

                TOther gamm_f = 1;
                // 边界定义
                auto A = (*bc->_FieldA3)[faceId];
                auto B = (*bc->_FieldB3)[faceId];
                auto C = (*bc->_FieldGamma3)[faceId];

                auto K = C / dcf;
                auto c = -sign * gamm_f * K / (C * (B + K)) * s;
                auto cf = B * c;
                auto vf = A * c;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }

        return true;
    }
    // 边界为第三类边界条件计算模态时，只处理常数项，该项并入刚度矩阵
    template <class TOther>
    bool MixedModal(BC<TCase, TField> *bc, const FieldBase<TCase, TOther> *gamma)
    {
        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;

        auto flux = this->_Equation->_Flux;
        long from = bc->start;
        long to = bc->end;
        TField sign = 1;
        if (gamma != 0)
            for (auto faceId = from; faceId <= to; faceId++)
            {

                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];

                auto gamm_f = (*gamma->fFaceValues)[faceId];
                // 边界中的常数部分
                auto A=(bc->_Imp->A[0])[faceId];

                auto cf = sign * gamm_f *A * s;
                auto vf = 0.0;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }
        else
            for (auto faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];

                auto s = _Mesh->faceArea[faceId];

                TOther gamm_f = 1.0;
                // 边界中的常数部分
                auto A=(bc->_Imp->A[0])[faceId];

                auto cf = sign * gamm_f *A * s;
                auto vf = 0.0;
                auto tf = cf * (*_phiField->fValues)[fC] + vf;

                (*flux->FluxC_f)[faceId] += cf;
                (*flux->FluxV_f)[faceId] += vf;
                (*flux->FluxT_f)[faceId] += tf;
            }

        return true;
    }
};
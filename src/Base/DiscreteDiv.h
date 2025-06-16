#pragma once
using namespace std;
#include "PhysicsBC.h"

template <class TCase, class TField>
class DiscreteBase;

template <class TCase, class TField>
class FieldBase;

template <class TCase, class TField>
class DiscreteDiv : public DiscreteBase<TCase, TField>
{
private:
    const FieldBase<TCase, TField> *_phiField;

public:
    DiscreteDiv(EquationBase<TCase, TField> *EQ)
        : DiscreteBase<T>(EQ)
    {
    }
    ~DiscreteDiv() {}
    // ∇(Γϕ) :其中 Γ为送入的场,ϕ为待求解场phiField
    // 将其转化为 sigma_f[(S_f.Γ_f)*Φ_f]
    // sigma_f 对网格cell的所有面f进行求和
    //_f为该面上的值
    // 1:Γ为矢量场 ，ϕ 均为标量场
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const FieldBase<TCase, vector3D<TOther>> *gamma)
    {
        _phiField = phiField;

        auto _Case = this->_Equation->_Case;
        auto _Mesh = _Case->_Mesh;
        if (!this->_Equation->_Signal && gamma != 0)
        {
            // 计算Γ的界面值
            if (!this->template CalcuFaceValuesVec<TOther>(gamma))
                return false;

            // 下次不需要计算
            this->_Equation->_Signal = true;
        }
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
    // div ((rho*U)*hf) :其中 rho为送入的scalarField, U为送入的矢量场，hf为target
    bool Discrete(FieldBase<TCase, TField> *scalarField,
                  FieldBase<TCase, TField> *vectorField1,
                  FieldBase<TCase, TField> *vectorField2,
                  FieldBase<TCase, TField> *vectorField3)
    {
        auto psiField = this->_Equation->_Case->GetScalarFieldFace("mdot");
        auto _Mesh = this->_Equation->_Case->pMesh;
        // 计算交界面值
        if (!this->_Equation->_Signal)
        {
            InterpolateFaceValue interpolate(this->_Equation->_Case, FaceValueSchemeType::Central);
            // 密度
            interpolate._InterpolateInPlace(scalarField);
            // 速度场
            interpolate._InterpolateInPlace(vectorField1, vectorField2, vectorField3);

            // 通量
            CalcuFlux calcu(this->_Equation->_Case);
            calcu.Calculate(scalarField, vectorField1, vectorField2, vectorField3, psiField->fFaceValues);

            // 下次不需要计算
            this->_Equation->_Signal = true;
        }

        // 内部面
        Inner();

        // 边界面上的处理,
        for (size_t i = 0; i < _Mesh->FaceBD.size(); i++)
        {
            auto _bc = _Mesh->faceSections[_Mesh->FaceBD[i]];
            auto bc = this->_Equation->_Case->GetBC(_bc);
            if (bc == NULL)
                throw "边界条件未定义";
            else
                Boundary(_bc, bc);
        }
        // 其他格式需要修正
        if (this->_Equation->pEquationConfigue->dmsFaceValueConvection != DiscreteMethodSchemeType::Upwind)
        {
            Correct();
        }

        return true;
    }

    void Inner()
    {
        auto flux = this->_Equation->_Flux[0];
        auto psiField = this->_Equation->_Case->GetScalarFieldFace("mdot");
        auto _Mesh = this->_Equation->_Case->pMesh;
        auto scalarFieldMain = this->_Equation->targetField;
        T sign = 1;

        size_t nFace = _Mesh->inFaces.size();
        for (size_t i = 0; i < nFace; i++)
        // foreach (var faceId in _Mesh.inFaces)
        {
            auto faceId = _Mesh->inFaces[i];
            // auto f = _Mesh->faces[faceId];
            auto fC = _Mesh->nextCell[faceId];
            auto fF = _Mesh->prevCell[faceId];
            auto psi = psiField->fFaceValues[faceId];
            auto cf = sign * fmax(psi, 0);
            auto ff = -sign * fmax(-psi, 0);
            auto vf = 0;
            double tf = cf * scalarFieldMain->fValues[fC] + ff * scalarFieldMain->fValues[fF] + vf;

            flux->FluxC_f[faceId] += cf;
            flux->FluxF_f[faceId] += ff;
            flux->FluxV_f[faceId] += vf;
            flux->FluxT_f[faceId] += tf;
        }
    }
    bool Boundary(tuple<long, long, int, int> bdFace, PhysicsBC<TField> *bc)
    {
        switch (bc->_BCType)
        {
        case PhysicsBCType::Wall:
            return Wall(bdFace, bc);
        default:
            throw "该边界条件未实现";
        }
        return true;
    }
    bool Wall(tuple<long, long, int, int> bdFace, PhysicsBC<TField> *bc)
    {
        return true;
    }
    void Correct()
    {
        throw "correct未实现";
    }
};
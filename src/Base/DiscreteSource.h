#pragma once
using namespace std;

template <class TCase, class TField>
class DiscreteBase;

template <class TCase, class TField>
class FieldBase;

template <class TCase, class TField>
class DiscreteSource : public DiscreteBase<TCase, TField>
{
public:
    DiscreteSource(EquationBase<TCase, TField> *EQ)
        : DiscreteBase<TCase, TField>(EQ)
    {
    }
    ~DiscreteSource() {}

    // Γ∙ϕ:其中 Γ为送入的场,ϕ为待求解场
    //  将其转化为Γ_v∙ϕ_v*V_v
    //_f为该CELL上的值
    // 采用显性/隐性格式：节点p的系数矩阵为V_v ∙ max(Γ_v,0)
    // 源项为V_v ∙ min(Γ_v,0)
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const FieldBase<TCase, TOther> *gamma)
    {
        auto pMesh = this->_Equation->_Case->_Mesh;
        long n = pMesh->nCells;
        auto eqC = this->_Equation->_EquationConfigue;
        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (eqC.dmsSource == DiscreteMethodSchemeType::Explicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 对源项直接乘以体积
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * (*gamma->fValues)[i] * (*phiField->fLast)[i];
            }
        }
        else if (eqC.dmsSource == DiscreteMethodSchemeType::Implicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 全隐式
                (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * (*gamma->fValues)[i];
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
               (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * max((*gamma->fValues)[i], 0.0);
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * min((*gamma->fValues)[i], 0.0) * (*phiField->fLast)[i];
            }
        }

        return true;
    }
    // Γ∙ϕ:其中 Γ为送入的变量,ϕ为待求解场
    //  将其转化为Γ_v∙ϕ_v*V_v
    //_f为该CELL上的值
    // 采用显性/隐性格式：节点p的系数矩阵为V_v ∙ max(Γ_v,0)
    // 源项为V_v ∙ min(Γ_v,0)
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const TOther gamma)
    {
        auto pMesh = this->_Equation->_Case->_Mesh;
        long n = pMesh->nCells;
        auto eqC = this->_Equation->_EquationConfigue;
        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (eqC.dmsSource == DiscreteMethodSchemeType::Explicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 对源项直接乘以体积
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * gamma * (*phiField->fLast)[i];
            }
        }
        else if (eqC.dmsSource == DiscreteMethodSchemeType::Implicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 全隐式
                (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * gamma;
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * max(gamma, 0.0);
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * min(gamma, 0.0) * (*phiField->fLast)[i];
            }
        }

        return true;
    }
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TField> *phiField, const TOther *gamma)
    {
        auto pMesh = this->_Equation->_Case->_Mesh;
        long n = pMesh->nCells;
        auto eqC = this->_Equation->_EquationConfigue;
        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (eqC.dmsSource == DiscreteMethodSchemeType::Explicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 对源项直接乘以体积
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * gamma[i] * (*phiField->fLast)[i];
            }
        }
        else if (eqC.dmsSource == DiscreteMethodSchemeType::Implicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 全隐式
                (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * gamma[i];
            }
        }
        else
        {
            for (int i = 0; i < n; i++)
            {
                (*flux->FluxC_e)[i] += sign * pMesh->cellVolume[i] * max(gamma[i], 0.0);
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * min(gamma[i], 0.0) * (*phiField->fLast)[i];
            }
        }

        return true;
    }
    
    // 直接对数组操作
    template <class TOther>
    bool Discrete(TOther *gamma)
    {
        auto pMesh = this->_Equation->_Case->_Mesh;
        long n = pMesh->nCells;
        auto eqC = this->_Equation->_EquationConfigue;
        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (eqC.dmsSource == DiscreteMethodSchemeType::Explicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 对源项直接乘以体积
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * gamma[i];
            }
        }
        return true;
    }
    template <class TOther>
    bool Discrete(const FieldBase<TCase, TOther> *gamma)
    {
        auto pMesh = this->_Equation->_Case->pMesh;
        long n = pMesh->nCells;
        auto eqC = this->_Equation->_EquationConfigue;
        auto flux = this->_Equation->_Flux;
        TField sign = 1;
        if (eqC.dmsSource == DiscreteMethodSchemeType::Explicit)
        {
            for (int i = 0; i < n; i++)
            {
                // 对源项直接乘以体积
                (*flux->FluxT_e)[i] += sign * pMesh->cellVolume[i] * (*gamma->fValues)[i];
            }
        }
        return true;
    }
};
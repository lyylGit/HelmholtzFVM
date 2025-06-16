#pragma once

#include "EquationBase.h"

// #define _USE_MATH_DEFINES

// #include "math.h"
using namespace std;

template <class TCase>
class CaseHelmholtz;


// 定义液体密度
template <class TCase, class TField>
class EquationHelmholtz : public EquationBase<TCase, TField>
{
public:
    // ∇∙∇(p)+k^2p=-iωρ0q
    EquationHelmholtz(CaseBase<TCase> *pCase)
        : EquationBase<TCase, TField>(pCase)
    {
    }
    ~EquationHelmholtz() {}

    bool Discrete(int iType)
    {
        this->_Flux->Clear();

        auto pi = PI;
        // 目标场为复数压力场
        auto pCase = (CaseHelmholtz<TCase> *)(this->_Case);

        TCase c0 = pCase->c0;    // 声速
        auto freq = pCase->freq; // 计算频率

        auto P = pCase->getF_P();     // 压力
        auto Rho = pCase->getF_Rho(); // 密度
        auto Q = pCase->getF_Q();     // 源项

        auto omega = freq * 2.0 * pi;
        auto k0 = 2 * pi * freq / c0;
        // ∇∙∇(p)
        this->template Laplacian<TCase>(P);
        // k^2p
        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Implicit;
        this->template Source<TCase>(P, k0 * k0);

        // iωρ0q
        MemSegment<TCase, TField> iq(this->_Case, "", DataStoreOn::onCell);
        // i
        iq.SetValues(complex<TCase>(0, 1));
        // iω
        iq.MultiInplace(omega);
        // ρ0
        iq.MultiInplaceF(Rho);
        // q
        iq.MultiInplaceF(Q);
        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Explicit;
        this->template Source<TField>(iq.data.get());

        // 主变量
        this->targetField = P;

        // 对于混合边界条件，
        // auto bcLists = pCase->template GetBoundaryConditions<TField>(P);
        // for (auto _bc : bcLists)
        // {
        //     auto bc = (BC<TCase, TField> *)_bc.get();
        //     if (bc->_BCType == BoundaryConditionType::Mixed)
        //     {

        //     }
        // }

        return true;
    }
};
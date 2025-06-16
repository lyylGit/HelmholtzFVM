#pragma once

#include "EquationBase.h"

// #define _USE_MATH_DEFINES
#include "EigenFEAST.h"

// #include "math.h"
using namespace std;

template <class TCase>
class CaseHelmholtzConfigDef;

// 定义液体密度
template <class TCase, class TField>
class EquationHelmholtzConfigDef : public EquationBase<TCase, TField>
{
public:
    // ∇∙∇(p)+k^2p=-iωρ0q
    EquationHelmholtzConfigDef(CaseBase<TCase> *pCase)
        : EquationBase<TCase, TField>(pCase)
    {
    }
    ~EquationHelmholtzConfigDef() {}

    void Init()
    {
    }

    bool Discrete(int iType)
    {
        if (this->_Case->_ProjectConfig->caseType == defCaseType::HelmModal)
        {
            return calcuModal();
        }
        else if (this->_Case->_ProjectConfig->caseType == defCaseType::HelmFreq)
        {
            return calcuFreq();
        }

        return true;
    }
    // CFD
    bool calcuFreq()
    {
        this->_Flux->Clear();

        // 目标场为复数压力场
        auto pCase = (CaseHelmholtzConfigDef<TCase> *)(this->_Case);

        // 根据定义置场
        MemSegment<TCase, TCase> p0(pCase, "", DataStoreOn::onCell);
        pCase->InitField(p0, pCase->_ProjectConfig->p0);

        MemSegment<TCase, TCase> gamma_0(pCase, "", DataStoreOn::onCell);
        pCase->InitField(gamma_0, pCase->_ProjectConfig->gamma);

        MemSegment<TCase, TCase> gamma_p1(gamma_0);
        gamma_p1.MultiInplaceM(p0);
        gamma_p1.DividedInplace(1.0);

        auto P = pCase->getF_P();            // 压力
        auto Rho = pCase->getF_Rho();        // 密度的倒数
        auto heatSource = pCase->heatSource; // 热源源项
        auto Q = pCase->getF_Q();            // 源项
        auto freq = pCase->freqCurrent;      // 计算频率
        auto omg = freq * 2.0 * PI;

        // ∇∙(1/rho*∇(p))
        this->_EquationConfigue.dmsFaceValue = DiscreteMethodSchemeType::Central;
        this->template Laplacian<TCase>(P, Rho);

        // ω^2/(γ*p0) p
        MemSegment<TCase, TCase> gamma_omg2(gamma_p1);
        gamma_omg2.MultiInplace(omg * omg);

        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Implicit;
        this->template Source<TCase>(P, gamma_omg2.data.get());

        // iω (γ-1)/(γ*p0) q
        MemSegment<TCase, TField> iq(gamma_0);
        iq.AddInplace(-1.0);
        iq.MultiInplaceM(gamma_p1);
        // 释热源项
        if (heatSource != NULL)
        {
            MemSegment<TCase, TField> hs(pCase, "", DataStoreOn::onCell);
            // 首先计算G*e^(jφ)
            heatSource->CalcuFTFSourceValues(freq, heatSource.get(), &hs);
            // 然后乘以(γ-1)/(γ*p0)
            hs.MultiInplaceM(iq);

            long refId = pCase->refId;
            // 修改CP中的项
            heatSource->GradSource(P, pCase->nref, refId, &hs, this);
        }
        // iω
        TField iomg = 1i;
        iomg *= omg;
        iq.MultiInplace(iomg);
        // q
        iq.MultiInplaceF(Q);
        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Explicit;
        this->template Source<TField>(iq.data.get());

        // 主变量
        this->targetField = P;
        return true;
    }
    // 模态
    bool calcuModal()
    {
        this->_Flux->Clear();

        // 目标场为复数压力场
        auto pCase = (CaseHelmholtzConfigDef<TCase> *)(this->_Case);

        // 根据定义置场
        MemSegment<TCase, TCase> p0(pCase, "", DataStoreOn::onCell);
        pCase->InitField(p0, pCase->_ProjectConfig->p0);

        MemSegment<TCase, TCase> gamma_0(pCase, "", DataStoreOn::onCell);
        pCase->InitField(gamma_0, pCase->_ProjectConfig->gamma);

        gamma_0.MultiInplaceM(p0);
        gamma_0.DividedInplace(1.0);

        auto P = pCase->getF_P();            // 压力
        auto Rho = pCase->getF_Rho();        // 密度的倒数
        auto heatSource = pCase->heatSource; // 热源源项

        // ∇∙∇(p)
        // 如果有与k无关的声阻抗系数，需要并入Stiffnes Matrix
        this->_Flux->Clear();
        this->_EquationConfigue.dmsFaceValue = DiscreteMethodSchemeType::Central;
        this->template Laplacian<TCase>(P, Rho);
        // Stiffnes Matrix
        map<long, map<long, TField>> K;
        auto dummy = make_shared<TField[]>(pCase->_Mesh->nCells);
        this->_Flux->AssembleMatrix(K, dummy.get(), dummy.get());

        // k^2p
        this->_Flux->Clear();
        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Implicit;
        this->template Source<TCase>(P, gamma_0.data.get()); // 1.0 / (gamma_0 * p_0));
        // Mass Matrix
        map<long, map<long, TField>> M;
        this->_Flux->AssembleMatrix(M, dummy.get(), dummy.get());

        // 组边界上的除了常数项的部分。
        map<int, map<long, map<long, TField>>> CP;
        CP[0] = K; // 加入 Stiffnes Matrix

        Squeeze(M);
        CP[2] = M; // 加入Mass Matrix

        map<long, map<long, TField>> p1;
        for (long i = 0; i < pCase->_Mesh->nCells; i++)
            p1[i][i] = 0.0;
        CP[1] = p1;

        // 频率变化的阻抗边界，注意此处为角频率拟合
        int maxK = 3;
        for (auto _bcc : pCase->_BCsNeedCalcu)
        {
            BC<TCase, complex<TCase>> *_bc = (BC<TCase, complex<TCase>> *)(_bcc.get());
            //
            maxK = max(maxK, _bc->_Imp->N);
            for (int i = 1; i < _bc->_Imp->N; i++)
            {
                auto A = &(_bc->_Imp->A[i]);
                for (long faceId = _bc->start; faceId <= _bc->end; faceId++)
                {
                    auto fC = pCase->_Mesh->ownerCell[faceId];
                    auto s = pCase->_Mesh->faceArea[faceId];
                    auto gamm_f = (*Rho->fFaceValues)[faceId]; // 密度的倒数 1/rho
                    auto cf = gamm_f * (*A)[faceId] * s;
                    CP[i][fC][fC] += cf;
                }
            }
        }

        // 频率变化的源项
        if (heatSource != NULL)
        {
            maxK = max(maxK, heatSource->N);
            long refId = pCase->refId;
            // 修改CP中的项
            heatSource->GradSource(P, pCase->nref, refId, CP);
        }

        int fpm[64];
        double center = pCase->_ProjectConfig->helmModal.center * 2 * PI; // 0-600Hz
        double Emid[2] = {center, 0.0}, r = pCase->_ProjectConfig->helmModal.r * 2 * PI;
        int M0 = pCase->_ProjectConfig->helmModal.M0;
        EigenParam<TCase> param(pCase->_Mesh->nCells, M0, maxK, fpm, Emid, r, pCase->_Mesh);
        param.fpm[0] = pCase->_ProjectConfig->helmModal.output;
        param.fpm[2] = pCase->_ProjectConfig->helmModal.err; // 1e-4;
        param.fpm[5] = pCase->_ProjectConfig->helmModal.criteria;
        param.fpm[42] = pCase->_ProjectConfig->helmModal.it; // 迭代求解器=1
        param.fpm[15] = 0;
        
        ::_GeneralEigenZ<TCase, TField>(CP, param);

        param.printInfo();
        if (param.info == 0)
        {
            string savepath = pCase->_ProjectConfig->outPath + "/"+pCase->_ProjectConfig->caseName+"_";
            param.printEigenVector(savepath, -1, false);
        }
        // 主变量
        this->targetField = P;

        return true;
    }
};
#pragma once

#include "EquationBase.h"

// #define _USE_MATH_DEFINES
#include "EigenFEAST.h"

// #include "math.h"
using namespace std;

template <class TCase>
class CaseHelmholtzModal;

// 定义液体密度
template <class TCase, class TField>
class EquationHelmholtzModal : public EquationBase<TCase, TField>
{
public:
    // ∇∙∇(p)+k^2p=-iωρ0q
    EquationHelmholtzModal(CaseBase<TCase> *pCase)
        : EquationBase<TCase, TField>(pCase)
    {
    }
    ~EquationHelmholtzModal() {}

    bool Discrete(int iType)
    {
        this->_Flux->Clear();

        // 目标场为复数压力场
        auto pCase = (CaseHelmholtzModal<TCase> *)(this->_Case);

        TCase c0 = pCase->c0; // 声速

        auto P = pCase->getF_P();     // 压力
        auto Rho = pCase->getF_Rho(); // 密度
        auto Q = pCase->getF_Q();     // 源项

        // ∇∙∇(p)
        // 如果有与k无关的声阻抗系数，需要并入Stiffnes Matrix
        this->_Flux->Clear();
        this->template Laplacian<TCase>(P);
        // Stiffnes Matrix
        map<long, map<long, TField>> K;
        auto dummy = make_shared<TField[]>(pCase->_Mesh->nCells);
        this->_Flux->AssembleMatrix(K, dummy.get(), dummy.get());

        // k^2p
        this->_Flux->Clear();
        this->_EquationConfigue.dmsSource = DiscreteMethodSchemeType::Implicit;
        this->template Source<TCase>(P, 1.0);
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
                    auto gamm_f = 1.0; // 没有gamma
                    auto cf = gamm_f * (*A)[faceId] * s;
                    CP[i][fC][fC] += cf;
                }
            }
        }

        int fpm[64];
        double Emid[2] = {sqrt(50.0), 0.0}, r = sqrt(50.0);
        int M0 = 180;
        EigenParam<TCase> param(pCase->_Mesh->nCells, M0, maxK, fpm, Emid, r, pCase->_Mesh);
        param.fpm[0] = 1;
        param.fpm[2] = 6; // 1e-6;
        // param.fpm[42] = 1; // 迭代求解器

        ::_GeneralEigenZ<TCase, TField>(CP, param);

        param.printInfo();
        param.printEigenVector("d:\\", 10, true);
        param.printEigenMatrix("d:\\modalR.txt");
        param.printEigenMatrix("d:\\modalL.txt",false);
        // param.printEigenMatrixBin("d:\\modal.bin", false);

        // ::WritteArray(K,pCase->_Mesh->nCells,"A0.txt");
        // ::WritteArray(p1,pCase->_Mesh->nCells,"A1.txt");
        // ::WritteArray(M,pCase->_Mesh->nCells,"A2.txt");
        // 主变量
        this->targetField = P;

        return true;
    }
};
#pragma once

#include <memory>
#include "MemSegment.h"

using namespace std;
// 单位体积激励源定义，采用FTF进行定义
// q=A0+A1ω+A2ω^2+A3ω^3+...
// ω为角频率：ω=2 * PI * freq
template <class TCase, class TField>
class FTFSource
{
public:
    int N; // 阶数+1
    CaseBase<TCase> *_Case;
    // 计算返回的A0,A1,A2,A3....
    shared_ptr<MemSegment<TCase, TField>[]> A;

    // 计算具体的释热脉动值，通常与ω有关
    function<void(double, FTFSource *, void *)> CalcuFTFSourceValues;
    // 计算系数，用来计算特征值问题
    function<void(FTFSource *)> CalcuFTFSourceCoeffs;

public:
    FTFSource(int _N, CaseBase<TCase> *pCase) : N(_N), _Case(pCase)
    {
        if (pCase->_ProjectConfig->caseType == defCaseType::HelmModal)
        {
            A = make_shared<MemSegment<TCase, TField>[]>(N);
            for (int i = 0; i < _N; i++)
                A[i].SetInfo(pCase, DataStoreOn::onCell);
        }
    }
    // 计算释热源项中的 grad（pref）* nref
    // Γ∙∇ϕ
    //  将其转化为 sigma_f[(Γ_f∙S_f*ϕ_f)]
    //  sigma_f 对网格cell的所有面f进行求和
    //_f为该面上的值，其中S_f为面积矢量
    // 注意，该描述中，源项在方程右侧，因此要加负号
    bool GradSource(const FieldBase<TCase, TField> *phiField, const vector3D<TCase> &gamma_f, long refId, map<int, map<long, map<long, TField>>> &CP)
    {
        // 参考单元的面
        map<long, TField> coeff;
        auto _Mesh = _Case->_Mesh;
        auto faceLists = _Mesh->_Cells[refId];
        for (auto faceId : faceLists)
        {
            auto fF = _Mesh->neighborCell[faceId];
            auto fC = _Mesh->ownerCell[faceId];
            // 采用中心差分
            auto gc = _Mesh->gc[faceId];

            auto S = _Mesh->surfaceVector[faceId];

            auto cf = gamma_f.dot(S);
            // 参考单元的面全部是内部面
            if (fF >= 0)
            {
                auto ff = -cf * (1 - gc);
                coeff[fF] += ff;
            }
            cf = -cf * gc;
            // if (fC == refId)
            coeff[fC] += cf;
        }
        // 对所有的cell 修正系数，不考虑释热区对参考区的影响
        for (int i = 0; i < N; i++)
        {
            for (long j = 0; j < _Mesh->nCells; j++)
            {
                if (abs(A[i][j]) > 0)
                {
                    for (auto it = coeff.begin(); it != coeff.end(); it++)
                    {
                        CP[i][j][it->first] += it->second * A[i][j];
                        // 其他影响
                        //  CP[i][it->first][j] -= it->second*A[i][j];
                    }
                }
            }
        }

        return true;
    }
    template <class TOther>
    bool GradSource(const FieldBase<TCase, TField> *phiField, const vector3D<TCase> &normal_f, long refId, MemSegment<TCase, TOther> *gamma, EquationBase<TCase, TField> *pEq)
    {
        // 参考单元的面
        map<long, TField> coeff;
        auto _Mesh = _Case->_Mesh;
        auto faceLists = _Mesh->_Cells[refId];
        for (auto faceId : faceLists)
        {
            auto fF = _Mesh->neighborCell[faceId];
            auto fC = _Mesh->ownerCell[faceId];
            // 采用中心差分
            auto gc = _Mesh->gc[faceId];

            auto S = _Mesh->surfaceVector[faceId];

            auto cf = normal_f.dot(S);
            // 参考单元的面全部是内部面
            if (fF >= 0)
            {
                // 负号是因为源项在方程右侧，移到左侧产生负号，下面相同
                auto ff = -cf * (1 - gc);
                coeff[fF] += ff;
            }
            cf = -cf * gc;
            coeff[fC] += cf;
        }

        auto flux = pEq->_Flux;
        // 对所有的cell 修正系数，不考虑释热区对参考区的影响

        for (long j = 0; j < _Mesh->nCells; j++)
        {
            if (abs((*gamma)[j]) > 0)
            {
                for (auto it = coeff.begin(); it != coeff.end(); it++)
                {
                    flux->SourceC[j][it->first] += it->second * (*gamma)[j];
                    flux->IndexC.push_back({j, it->first});
                }
            }
        }

        return true;
    }
};
#pragma once

#include "MemSegment.h"

template <class T>
class CaseBase;
template <class T>
class Mesh3;

template <class TCase, class TField>
class Flux //: public staticFunctions
{

public:
    CaseBase<TCase> *_Case;
    // 网格单元个数,交界面个数//,顶点数:结冰计算用
    long nCells = 0, nFaces = 0; //, nVertex = 0;
    // 网格
    Mesh3<TCase> *_Mesh;
    // Cell[] _Cells;

    // P109
    // 面通量: C_f*phi_C+F_f*phi_F+V_f
    // FluxC_f:Cell_C在f面上的通量系数
    // FluxF_f:Cell_F在f面上的通量系数
    // FluxV_f:非线性部分
    // FluxT_f:面的总通量
    shared_ptr<MemSegment<TCase, TField>> FluxC_f, FluxF_f, FluxV_f, FluxT_f;
    // 面上的通量值
    shared_ptr<MemSegment<TCase, TField>> mdot_f, mdot_f_prev;
    // 控制体通量,源项部分 C_e*phi_C+V_e
    //  FluxC_e:Cell_C的体积分系数
    //  FluxV_e:非线性部分
    //  FluxT_e 体的总通量
    shared_ptr<MemSegment<TCase, TField>> FluxC_e, FluxC_e_old, FluxV_e, FluxT_e;

    // 系数矩阵
    // unordered_map<long, unordered_map<long, TField>> MatrixA;
    map<long, map<long, TField>> MatrixA;
    vector<pair<long, long>> IndexA;

    // 针对直接采用网格单元间关系构建的附加信息，例如FTF产生的对参考点的依赖
    map<long, map<long, TField>> SourceC;
    vector<pair<long, long>> IndexC;

    shared_ptr<TField[]> VectorB, AC_old; // MatrixA的前一步对角线分量
public:
    Flux(CaseBase<TCase> *pCase)
    {
        _Case = pCase;
        _Mesh = pCase->_Mesh;
        nCells = _Mesh->nCells;
        nFaces = _Mesh->nFaces;

        Init();
    }
    ~Flux() {}

private:
    void Init()
    {

        VectorB = make_shared<TField[]>(nCells);
        AC_old = make_shared<TField[]>(nCells);

        // 按面
        FluxC_f = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onFace);
        FluxF_f = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onFace);
        FluxV_f = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onFace);
        FluxT_f = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onFace);

        // 按cell
        FluxC_e = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onCell);
        FluxC_e_old = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onCell);
        FluxV_e = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onCell);
        FluxT_e = make_shared<MemSegment<TCase, TField>>(_Case, DataStoreOn::onCell);
    }

public:
    void Clear()
    {
        Zero(MatrixA, IndexA);
        Zero(SourceC, IndexC);
        // 清除内存
        memset(VectorB.get(), 0, nCells * sizeof(TField));
        memset(AC_old.get(), 0, nCells * sizeof(TField));

        FluxC_f->Zero();
        FluxF_f->Zero();
        FluxV_f->Zero();
        FluxT_f->Zero();
        FluxC_e->Zero();
        FluxC_e_old->Zero();
        FluxV_e->Zero();
        FluxT_e->Zero();
    }
    void AssembleMatrix()
    {
        AssembleMatrix(MatrixA, VectorB.get(), AC_old.get());
    }
    void AssembleMatrix(map<long, map<long, TField>> &_MatrixA, TField *_VectorB, TField *_AC_old)
    {
        // 内部面
        size_t n = _Mesh->inFaces.size();
        for (size_t i = 0; i < n; i++)
        {
            long faceId = _Mesh->inFaces[i];
            long C = _Mesh->ownerCell[faceId],
                 F = _Mesh->neighborCell[faceId];

            _MatrixA[C][C] += (*FluxC_f)[faceId];
            _MatrixA[C][F] += (*FluxF_f)[faceId];
            VectorB[C] += -(*FluxT_f)[faceId];

            _MatrixA[F][F] += -(*FluxF_f)[faceId];
            _MatrixA[F][C] += -(*FluxC_f)[faceId];
            _VectorB[F] += (*FluxT_f)[faceId];
        }

        // 边界面
        n = _Mesh->bdFaces.size();
        for (size_t i = 0; i < n; i++)
        {
            long faceId = _Mesh->bdFaces[i];
            int C = _Mesh->ownerCell[faceId];

            _MatrixA[C][C] += (*FluxC_f)[faceId];
            _VectorB[C] += -(*FluxT_f)[faceId];
        }

        for (long C = 0; C < nCells; C++)
        {
            _MatrixA[C][C] += (*FluxC_e)[C];
            _VectorB[C] += -(*FluxT_e)[C];

            _AC_old[C] += (*FluxC_e_old)[C];
        }
        // 源项修正
        for (auto ij : IndexC)
        {
            auto i = get<0>(ij);
            auto j = get<1>(ij);
            _MatrixA[i][j] += SourceC[i][j];
        }
    }
    // 求解 A x= B
    bool IterateEquation(TField *Result)
    {
        //
        IndexMatrix(MatrixA, IndexA);
        // 输出矩阵
        // _OutputMatrix(MatrixA, VectorB.get(), nCells);

        bool b = true;
        b = _IterateEquation<TField>(MatrixA, IndexA, VectorB.get(), Result, nCells, 16);

        return b;
    }
    // 求解 CT*A*C x= CT*B,C为N*M矩阵，M<N
    bool IterateEquation(TField *Result, TField *C, long M)
    {
        //
        IndexMatrix(MatrixA, IndexA);

        bool b = true;
        b = _IterateEquation<TField>(MatrixA, IndexA, VectorB.get(), Result, nCells, C, M, 16);

        return b;
    }
};

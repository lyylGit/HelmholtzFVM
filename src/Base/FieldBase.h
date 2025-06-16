#pragma once

#include "MemSegment.h"
#include <functional>
#include "staticFunctions.h"
#include "Mesh3.h"
#include "Types.h"
#include "CaseBase.h"
#include "vector3d.h"

using namespace std;

template <class TCase, class TField>
class FieldBase
{
public:
    // 脚本描述中所使用的符号
    string Name;
    // 类型:标量或矢量
    FieldType fType;
    // 类型:物性或待求解变量
    VariableType vType;

    // 体积值:当前值,前一时刻,前前时刻值
    shared_ptr<MemSegment<TCase, TField>> fValues, fValuesPrev, fValuesPrevPrev;
    // 交界面值:当前值,前一时刻,前前时刻值
    shared_ptr<MemSegment<TCase, TField>> fFaceValues, fFaceValuesPrev, fFaceValuesPrevPrev;
    // 上一步迭代值
    shared_ptr<MemSegment<TCase, TField>> fLast;

    // 关联的Case
    CaseBase<TCase> *_Case;

public:
    FieldBase(CaseBase<TCase> *pCase)
    {
        _Case = pCase;
    }
    FieldBase(CaseBase<TCase> *pCase, string name)
    {
        _Case = pCase;
        Name = name;
    }
    ~FieldBase()
    {
    }

public:
    // 初始化，设置数组空间
    void SetDimension()
    {
        long ncells = _Case->_Mesh->nCells,
             nfaces = _Case->_Mesh->nFaces;
        bool stability = _Case->Stability;

        fValues = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onCell);
        fFaceValues = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onFace);
        fLast = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onCell);
        if (!stability)
        {
            fValuesPrev = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onCell);
            fValuesPrevPrev = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onCell);

            fFaceValuesPrev = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onFace);
            fFaceValuesPrevPrev = make_shared<MemSegment<TCase, TField>>(_Case, Name, DataStoreOn::onFace);
        }
    }
    // 从Case中查找其他场
    FieldBase *GetField(string name)
    {
        // 没检查,自己保证不会为null
        if (_Case != NULL)
            return _Case->fields[name];
        else
            return NULL;
    }
    
    // 设置全场为某值
    void SetCellValue(TField f, bool bIni)
    {
        // cell的值
        fValues->SetValues(f);
        fLast->SetValues(f);
        if (bIni && !_Case->Stability)
            fValuesPrev->SetValues(f);

        // face的值
        fFaceValues->SetValues(f);

        if (bIni && !_Case->Stability)
            fFaceValuesPrev->SetValues(f);
    }
     void SetCellValue(TField* f, bool bIni)
    {
        // cell的值
        fValues->SetValues(f);
        fLast->SetValues(f);
        if (bIni && !_Case->Stability)
            fValuesPrev->SetValues(f);

        // face的值,需要差值
        CalcuFaceValue(bIni);
      
    }
    template <class TOther>
    void SetCellValue(MemSegment<TCase, TOther>& f, bool bIni)
    {
        // cell的值
        fValues->SetValues(f);
        fLast->SetValues(f);
        if (bIni && !_Case->Stability)
            fValuesPrev->SetValues(f);

        // face的值,需要差值
        CalcuFaceValue(bIni);
      
    }
    // 置某些单元的值，一般是边界
    void SetCellValue(TField f, vector<long> &cellIDs, bool bPrev)
    {
        if (cellIDs.empty())
            return;

        fValues->SetValues(f, cellIDs);
        fLast->SetValues(f, cellIDs);

        if (bPrev)
            fValuesPrev->SetValues(f, cellIDs);
    }
    // // 计算面上的值
    void CalcuFaceValue(bool bIni)
    {
        ::CalcuFaceValues(_Case->_Mesh, NULL, this);
        if (bIni)
            SetBCValueFromCellValue(*this);

        // 赋前值
        if (bIni && !_Case->Stability)
            fFaceValuesPrev->CopyData(*(fFaceValues.get()));
    }
    // 从场的值设置边界上的值
    void SetBCValueFromCellValue(FieldBase &f)
    {
        for (auto bd : _Case->_Mesh->FaceBD)
        {
            auto item = _Case->_Mesh->faceSections[bd];
            int istart = get<0>(item), iend = get<1>(item); //, bcType = get<2>(item);
            for (int faceId = istart; faceId <= iend; faceId++)
            {
                auto ownerCell = _Case->_Mesh->ownerCell[faceId];
                f.fFaceValues->data[faceId] = f.fValues->data[ownerCell];
            }
        }
    }

public:
    void Calculate(int time);

    // 由各个子类继承实现
    // bIni 是否是初始化
    virtual void SetCellValueOnCaseState(bool bIni){}
    virtual void UpdateCellValue(){}

public:
    // 采用中心格式计算交界面的值

public:
    
    // 算子：计算本Field在网格上的梯度值
    
    shared_ptr<vector3D<TField>[]>
    CalcuVolumeGrad()
    {
       
    }
};
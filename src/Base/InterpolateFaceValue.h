#pragma once
#include "Types.h"

template <class TCase,class TField>
class EquationBase;
template <class TCase,class TField>
class FieldBase;
template <class T>
class CaseBase;
template <class TCase,class TField>
class InterpolateFaceValue
{
private:
    /* data */
public:
    InterpolateFaceValue(CaseBase<TCase> *pCase,FaceValueSchemeType AvgType)
    {
    _Case = pCase;
    _AvgType = AvgType;
    }
    ~InterpolateFaceValue(){}

    CaseBase<TCase> *_Case;
    FaceValueSchemeType _AvgType;
   
   public:
   //标量的面插值
   void _InterpolateInPlace(FieldBase<TCase,TField>* _Field)
   {
     auto mesh = _Case->pMesh;
    if (_AvgType == FaceValueSchemeType::Central)
    {
        size_t n = mesh->inFaces.size();
        for (size_t i = 0; i < n; i++)
        {
            long faceId = mesh->inFaces[i];
            auto fC = mesh->nextCell[faceId];
            auto fF = mesh->prevCell[faceId];
            auto gc = mesh->gc[faceId];
            _Field->fFaceValues[faceId] = gc * _Field->fValues[fC] + (1 - gc) * _Field->fValues[fF];
        }
    }
   }
   //矢量的面插值
   void _InterpolateInPlace(FieldBase<TCase,TField>* _vField1,FieldBase<TCase,TField>* _vField2,FieldBase<TCase,TField>* _vField3)
   {
     auto mesh = _Case->pMesh;
    if (_AvgType == FaceValueSchemeType::Central)
    {
        size_t n = mesh->inFaces.size();
        for (size_t i = 0; i < n; i++)
        {
            long faceId = mesh->inFaces[i];
            auto fC = mesh->nextCell[faceId];
            auto fF = mesh->prevCell[faceId];
            auto gc = mesh->gc[faceId];
            _vField1->fFaceValues[faceId] = gc * _vField1->fValues[fC] + (1 - gc) * _vField1->fValues[fF];
            _vField2->fFaceValues[faceId] = gc * _vField2->fValues[fC] + (1 - gc) * _vField2->fValues[fF];
            _vField3->fFaceValues[faceId] = gc * _vField3->fValues[fC] + (1 - gc) * _vField3->fValues[fF];
        }
    }
   }
};
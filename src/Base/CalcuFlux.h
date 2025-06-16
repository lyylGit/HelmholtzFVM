#pragma once
using namespace std;
#include "PhysicsBC.h"

template <class TCase, class TField>
class DiscreteBase;

template <class TCase, class TField>
class FieldBase;
template <class TCase, class TField>
class CalcuFlux
{
public:
    CaseBase<TCase> *_Case;

    CalcuFlux(CaseBase<TCase> *pCase)
    {
        _Case = pCase;
    }
    ~CalcuFlux() {};

public:
    bool Calculate(FieldBase<TCase, TField> *scalar, FieldBase<TCase, TField> *u1, FieldBase<TCase, TField> *u2, FieldBase<TCase, TField> *u3, TField *flux)
    {
        auto _Mesh = this->_Case->pMesh;

        Inner(scalar, u1, u2, u3, flux);
        // 边界面
        size_t n = _Mesh->FaceBD.size();
        for (size_t i = 0; i < n; i++)
        {
            auto bf = _Mesh->faceSections[_Mesh->FaceBD[i]];
            auto bc = this->_Case->GetBC(bf);
            Boundary(bf, bc, flux);
        }
        return true;
    }

    bool Inner(FieldBase<TCase, TField> *scalar, FieldBase<TCase, TField> *u1, FieldBase<TCase, TField> *u2, FieldBase<TCase, TField> *u3, TField *flux)
    {
         auto _Mesh = this->_Case->pMesh;
    size_t n = _Mesh->inFaces.size();
    for (size_t i = 0; i < n; i++)
    {
        auto faceId = _Mesh->inFaces[i];
        Vector fn(_Mesh->faceNormal, faceId);
        auto S = fn * _Mesh->faceArea[faceId];
        auto rho = scalar->fFaceValues[faceId];
        Vector U(u1->fFaceValues[faceId], u2->fFaceValues[faceId], u3->fFaceValues[faceId]);
        flux[faceId] = rho * U * S; // _U.fFaceValues[faceId].dot(S);
    }
    return true;
    }
    bool Boundary(tuple<long, long, int, int> bdFace, PhysicsBC<TField> *bc, TField *flux)
    {
         switch (bc->_BCType)
    {
    case PhysicsBCType::Wall:
        return Wall(bdFace, bc, flux);
    default:
        throw "该边界条件未实现";
    }
    }
    bool Wall(tuple<long, long, int, int> bdFace, PhysicsBC<TField> *bc, TField *flux)
    {
         long start = get<0>(bdFace), end = get<1>(bdFace);
    for (long faceId = start; faceId <= end; faceId++)
    {
        flux[faceId] = 0;
    }
    return true;
    }
};

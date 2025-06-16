#pragma once

using namespace std;

template <class TCase, class TField>
class CalcuVolumeGrad : public CalcuBase<TCase, TField>
{
    // 计算网格单元的梯度Grad
    // pp275：Green-Gauss Method
public:
    CalcuVolumeGrad(const FieldBase<TCase, TField> *pField, EquationConfigue *eqC)
        : CalcuBase<TCase, TField>(pField, eqC)
    {
    }

    ~CalcuVolumeGrad() {}

    // 返回值
    // 返回的是一个矢量场
    shared_ptr<vector3D<TField>[]> _Grad;

public:
    bool Calcu()
    {
        // 采用Gauss计算
        Gauss();

        return true;
    }

private:
    // 采用Gauss公式计算
    void Gauss()
    {
        auto _Field = this->_Field;
        auto _Case = _Field->_Case;
        auto _Mesh = _Case->_Mesh;

        long nCells = _Mesh->nCells;

        // 分配内存
        _Grad = make_shared<vector3D<TField>[]>(nCells);

        // Inner
        for (auto faceId : _Mesh->inFaces)
        {
            auto fF = _Mesh->neighborCell[faceId];
            auto fC = _Mesh->ownerCell[faceId];
            auto S = _Mesh->surfaceVector[faceId];
            auto gc = _Mesh->gc[faceId];
            auto phi_f = gc * (*_Field->fValues)[fC] + (1 - gc) * (*_Field->fValues)[fF];

            vector3D<TField> phi_f_S(phi_f * S.x, phi_f * S.y, phi_f * S.z);

            _Grad[fC] += phi_f_S;
            _Grad[fF] += -phi_f_S;
        }

        // Boundary
        for (auto _bcId : _Mesh->FaceBD)
        {
            BC<TCase, TField> *bc = (BC<TCase, TField> *)_Case->template GetBoundaryConditions<TField>(_Field, _bcId);
            if (bc == 0)
            {
                BC<TCase, TField> _bc("", _bcId);
                _bc.InitInfo(_Mesh->faceSections[_bcId]);
                _bc._BCType = BoundaryConditionType::UseCenter;
                BDGauss(&_bc);
            }
            else
                BDGauss(bc);
        }
        for (int iCell = 0; iCell < nCells; iCell++)
        {
            auto volumn = _Mesh->cellVolume[iCell];
            _Grad[iCell] /= volumn;
        }
    }
    // 计算边界对控制体梯度的影响
    void BDGauss(BC<TCase, TField> *bc)
    {
        auto _Field = this->_Field;
        auto _Case = _Field->_Case;
        auto _Mesh = _Case->_Mesh;

        // auto bdFace = bc->_BCFaces;
        long from = bc->start;
        long to = bc->end;
        switch (bc->_BCType)
        {
        // 边界上是值
        case BoundaryConditionType::Value:
        {
            for (int faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                // 直接采用面上的值
                auto f = (*_Field->fFaceValues)[faceId];
                vector3D<TField> grad(f * S.x, f * S.y, f * S.z);
                _Grad[fC] += grad;
            }
            break;
        }
        // 没有定义,界面上值采用中心值
        case BoundaryConditionType::UseCenter:
        {
            for (int faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                auto f = (*_Field->fValues)[fC];
                vector3D<TField> grad(f * S.x, f * S.y, f * S.z);
                _Grad[fC] += grad;
            }

            break;
        }
        // 界面上是梯度
        case BoundaryConditionType::Gradient:
        {
            for (int faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];
                auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                // 网格中心到边界面中心的距离
                auto d = _Mesh->dcf[faceId];
                // 插值得出边界通量,后面要乘以面积，所以此处先除以面积
                // auto uf = (*_Field->fFaceValues)[faceId] * d / s + (*_Field->fValues)[fC];
                auto uf = (*_Field->fFaceValues)[faceId] * d + (*_Field->fValues)[fC];
                _Grad[fC] += vector3D<TField>(uf * S.x, uf * S.y, uf * S.z);
            }
            break;
        }
        // 第三类边界条件
        case BoundaryConditionType::Mixed:
        {
            for (int faceId = from; faceId <= to; faceId++)
            {
                auto fC = _Mesh->ownerCell[faceId];
                // auto s = _Mesh->faceArea[faceId];
                auto S = _Mesh->surfaceVector[faceId];
                // 网格中心到边界面中心的距离
                auto dcf = _Mesh->dcf[faceId];
                // 边界定义
                auto A = (*bc->_FieldA3)[faceId];
                auto B = (*bc->_FieldB3)[faceId];
                auto C = (*bc->_FieldGamma3)[faceId];
                // 插值得出边界面上的值
                auto K = C / dcf;
                auto f = (K * (*_Field->fValues)[fC] - A) / (B + K);
                vector3D<TField> grad(f * S.x, f * S.y, f * S.z);
                _Grad[fC] += grad;
            }
            break;
        }
        }
    }
};
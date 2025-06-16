#pragma once
using namespace std;

template <class TCase,class TField>
class DiscreteBase;

template <class TCase,class TField>
class FieldBase;

template <class TCase,class TField>
class DiscreteDt : public DiscreteBase<TCase,TField>
{
public:
    DiscreteDt(EquationBase<TCase,TField> *EQ)
    : DiscreteBase<T>(EQ)
{
}
    ~DiscreteDt(){}

    bool Discrete(FieldBase<TCase,TField> *scalarField)
    {
         switch (this->_Equation->pEquationConfigue->dmsTransient)
    {
    case DiscreteSchemeTransient::BackwardEuler:
        return DiscreteBackEuler(scalarField);
    case DiscreteSchemeTransient::AdamsMoulton:
        return DiscreteAdamsMoulton(scalarField);
    default:
        this->_Equation->_Case->pMesh->ShowPrompt( "暂不支持");
        return false;
    }
    }

    bool DiscreteBackEuler(FieldBase<TCase,TField> *scalarField)
    {
        auto pMesh = this->_Equation->_Case->pMesh;
    //对系数矩阵进行赋值
    //var cells = _Case-pMesh;
    long n =pMesh->nCells;
    double dt = this->_Equation->pEquationConfigue->_DeltaT;


//针对矢量
    // if (_Target >= 0)
    // {
    //     var flux = _Equation.Flux[_Target];

    //     for (int iCell = 0; iCell < n; iCell++)
    //     {
    //         var ce = sign * scalarField.fValues[iCell] * cells[iCell].volume / dt;
    //         // var ce_old = -sign * scalarField.fFaceValuesPrev[iCell] * cells[iCell].volume / dt;
    //         var ce_old = -sign * scalarField.fValuesPrev[iCell] * cells[iCell].volume / dt;
    //         var te = /* ce * vectorFieldMain.fValues[iCell].x[_Target] +*/ ce_old * vectorFieldMain.fValuesPrev[iCell].x[_Target];

    //         flux.FluxC_e[iCell] += ce;
    //         flux.FluxC_e_old[iCell] += ce_old;
    //         flux.FluxT_e[iCell] += te;
    //     }
    // }
    // else
    // {
        auto flux =this->_Equation->_Flux[0];
        auto scalarFieldMain = this->_Equation->targetField;
        for (long iCell = 0; iCell < n; iCell++)
        {
            auto ce =  scalarField->fValues[iCell] *pMesh->cellVolume[iCell] / dt;
            auto ce_old = - scalarField->fValuesPrev[iCell] * pMesh->cellVolume[iCell]/ dt;
            auto te = /* ce * scalarFieldMain.fValues[iCell] +*/ ce_old * scalarFieldMain->fValuesPrev[iCell];
            flux->FluxC_e[iCell] += ce;
            flux->FluxC_e_old[iCell] += ce_old;
            flux->FluxT_e[iCell] += te;
        }
    // }

    return true;
    }
    bool DiscreteAdamsMoulton(FieldBase<TCase,TField> *scalarField)
    {
        auto pMesh = this->_Equation->_Case->pMesh;
    //对系数矩阵进行赋值
    //var cells = _Case-pMesh;
    long n =pMesh->nCells;
    double dt = this->_Equation->pEquationConfigue->_DeltaT;

    // if (_Target >= 0)
    // {
    //     var flux = _Equation.Flux[_Target];

    //     for (int iCell = 0; iCell < n; iCell++)
    //     {
    //         var ce = 1.5 * sign * scalarField.fValues[iCell] * cells[iCell].volume / dt;
    //         var ce_old = -2.0 * sign * scalarField.fFaceValuesPrev[iCell] * cells[iCell].volume / dt;
    //         var ve = 0.5 * sign * scalarField.fFaceValuesPrevPrev[iCell] * cells[iCell].volume * vectorFieldMain.fValuesPrevPrev[iCell].x[_Target] / dt;

    //         var te = ce * vectorFieldMain.fValues[iCell].x[_Target] + ce_old * vectorFieldMain.fValuesPrev[iCell].x[_Target] + ve;

    //         flux.FluxC_e[iCell] += ce;
    //         flux.FluxC_e_old[iCell] += ce_old;
    //         flux.FluxT_e[iCell] += te;
    //     }
    // }
    // else
    // {
        auto flux = this->_Equation->_Flux[0];
        auto scalarFieldMain = this->_Equation->targetField;

        for (long iCell = 0; iCell < n; iCell++)
        {
            auto ce = 1.5 *  scalarField->fValues[iCell] * pMesh->cellVolume[iCell] / dt;
            auto ce_old = -2.0 *  scalarField->fFaceValuesPrev[iCell] * pMesh->cellVolume[iCell] / dt;
            auto ve = 0.5 *  scalarField->fFaceValuesPrevPrev[iCell] * pMesh->cellVolume[iCell] * scalarFieldMain->fValuesPrevPrev[iCell] / dt;

            auto te = ce * scalarFieldMain->fValues[iCell] + ce_old * scalarFieldMain->fValuesPrev[iCell] + ve;

            flux->FluxC_e[iCell] += ce;
            flux->FluxC_e_old[iCell] += ce_old;
            flux->FluxT_e[iCell] += te;
        }
    // }
    return true;
    }
};
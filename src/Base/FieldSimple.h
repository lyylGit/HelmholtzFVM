#pragma once
#include "FieldBase.h"

using namespace std;

template <class TCase, class TField>
class FieldSimple : public FieldBase<TCase,TField>
{
public:
    FieldSimple(CaseBase<TCase> *_Case, string name, function<void()> fUpdate)
     : FieldBase<TCase,TField>(_Case)
    {
        this->Name = name;
        this->fType = FieldType::_Scalar;
        this->vType = VariableType::CAL;

        this->SetDimension();

        _fUpdate = fUpdate;
    }

    function<void()> _fUpdate;
    // 由各个子类继承实现
    // bIni 是否是初始化
    void SetCellValueOnCaseState(bool bIni){};
    void UpdateCellValue() { _fUpdate(); };

    
};
#pragma once

using namespace std;

template <class TCase, class TField>
class FieldBase;
template <class TCase, class TField>
class Flux;

template <class TCase, class TField>
class BC;

template <class TCase>
class CaseBase;

template <class TCase, class TField>
class EquationBase
{
public:
    // 方程名称
    string _Name;
    // 主要求解变量
    string _Variable;
    // 关联的Case
    CaseBase<TCase> *_Case;
    // 触发重新计算
    bool _Signal;

    EquationConfigue _EquationConfigue;

    // 通量集
    shared_ptr<Flux<TCase, TField>> _Flux;

    // 本方程面向的待求解场
    FieldBase<TCase, TField> *targetField;

public:
    // 构造器
    EquationBase(CaseBase<TCase> *pCase)
    {
        _Signal = false;
        _Case = pCase;
        _Flux = make_shared<Flux<TCase, TField>>(pCase);
    }
    ~EquationBase() {}

    virtual void Init() {};

public:
    // 离散
    // ∇∙(Γ∇ϕ)
    template <class TOther>
    bool Laplacian(FieldBase<TCase, TField> *phi, FieldBase<TCase, TOther> *gamma = 0)
    {
        DiscreteLaplacian<TCase, TField> laplacian(this);

        laplacian.template Discrete<TOther>(phi, gamma);

        return true;
    }

    // 离散源项
    // Γ∙ϕ
    template <class TOther>
    bool Source(FieldBase<TCase, TField> *phi, FieldBase<TCase, TOther> *gamma = 0)
    {
        DiscreteSource<TCase, TField> source(this);
        if (gamma != 0)
            source.template Discrete<TOther>(phi, gamma);
        else
            source.template Discrete<TOther>(phi, 1.0);

        return true;
    }
   
    // 离散源项
    // a∙ϕ
    template <class TOther>
    bool Source(FieldBase<TCase, TField> *phi, const TOther a)
    {
        DiscreteSource<TCase, TField> source(this);
        source.template Discrete<TOther>(phi, a);

        return true;
    }
    // A∙ϕ
    template <class TOther>
    bool Source(FieldBase<TCase, TField> *phi, const TOther *a)
    {
        DiscreteSource<TCase, TField> source(this);
        source.template Discrete<TOther>(phi, a);

        return true;
    }
    // A
    template <class TOther>
    bool Source(TOther *a)
    {
        DiscreteSource<TCase, TField> source(this);
        source.template Discrete<TOther>(a);
        return true;
    }
    // A
    template <class TOther>
    bool Source(FieldBase<TCase, TOther> *a)
    {
        DiscreteSource<TCase, TField> source(this);
        source.template Discrete<TOther>(a);
        return true;
    }

    virtual bool Discrete(int itype) { return true; }
    // 离散各项
    // 瞬态项:例如 rho*∂hf/∂t; 其中rho为scale送进去,hf为本方程的targetField
    // bool Dt(FieldBase<T> *scalar);
    // // 对流项:例如 div(rhow*u*hf)；其中 rohw为scalar,hf为本方程的targetField，u的三个分量送入
    // bool Div(FieldBase<T> *scalar,
    //          FieldBase<T> *vectorField1, FieldBase<T> *vectorField2, FieldBase<T> *vectorField3);
    // // 源项，目前支持显式格式
    // bool Source(FieldBase<T> *scalar);
    // bool Source(T *scalar);

public:
    // 迭代
    bool IterateEquation()
    {
        _Flux->AssembleMatrix();
        // TField *R = new TField[_Case->_Mesh->nCells];
        bool b = true;
        b = _Flux->IterateEquation(targetField->fValues->data.get());

        // delete[] R;

        return b;
    }

    // 迭代
    bool IterateEquation(TField *C, long M)
    {
        _Flux->AssembleMatrix();
        // TField *R = new TField[_Case->_Mesh->nCells];
        bool b = true;
        b = _Flux->IterateEquation(targetField->fValues->data.get(), C, M);

        // delete[] R;

        return b;
    }
};

#pragma once
#include <iostream>
#include <functional>
#include <cstring>
#include "Mesh3.h"
#include "Types.h"
#include "CaseBase.h"

using namespace std;

#define _NODISCARD [[nodiscard]]

template <class TCase, class TField>
class MemSegment
{
    // 本类主要封装场数据和操作
public:
    // 脚本描述中所使用的符号
    string Name;
    // 关联的Case
    CaseBase<TCase> *_Case;
    // 数据长度
    long length;
    // 存储于点、面、体上？
    DataStoreOn dataStoreOn;
    // 数据
    shared_ptr<TField[]> data;

    MemSegment()
    {
    }

    MemSegment(CaseBase<TCase> *pCase, DataStoreOn _dataStoreOn)
        : MemSegment(pCase, "", _dataStoreOn)
    {
    }
    MemSegment(CaseBase<TCase> *pCase, string name, DataStoreOn _dataStoreOn)
        : Name(name), _Case(pCase), dataStoreOn(_dataStoreOn), data(0)
    {
        SetInfo(pCase, _dataStoreOn);
    }
    MemSegment(CaseBase<TCase> *pCase, string name, string fromname, DataStoreOn _dataStoreOn)
        : Name(name), _Case(pCase), dataStoreOn(_dataStoreOn)
    {
        switch (_dataStoreOn)
        {
        case DataStoreOn::onFace:
            length = pCase->_Mesh->nFaces;
            data = pCase->_Mesh->faceData[fromname];
            break;
        case DataStoreOn::onCell:
            length = pCase->_Mesh->nCells;
            data = pCase->_Mesh->cellData[fromname];
            break;

        default:
            break;
        }
    }
    template <class TOther>
    MemSegment(const MemSegment<TCase, TOther> &clone)
    {
        DeepCopy(clone);
    }

    ~MemSegment()
    {
    }
    void SetInfo(CaseBase<TCase> *pCase, DataStoreOn _dataStoreOn)
    {
        _Case = pCase;
        dataStoreOn = _dataStoreOn;

        switch (_dataStoreOn)
        {
        case DataStoreOn::onFace:
            length = pCase->_Mesh->nFaces;
            break;
        case DataStoreOn::onCell:
            length = pCase->_Mesh->nCells;
            break;
        default:
            break;
        }
        data = make_shared<TField[]>(length);
        memset(data.get(), 0, length * sizeof(TField));
    }
    void Zero()
    {
        memset(data.get(), 0, length * sizeof(TField));
    }
    template <class TOther>
    void DeepCopy(const MemSegment<TCase, TOther> &clone)
    {
        Name = clone.Name;
        _Case = clone._Case;
        dataStoreOn = clone.dataStoreOn;
        length = clone.length;
        data = make_shared<TField[]>(length);
        // memcpy(data.get(), clone.data.get(), length * sizeof(TField));
        for (long i = 0; i < length; i++)
            data[i] = clone.data[i];

        // cout << "deepcopy\n";
    }
    void ShallowCopy(const MemSegment &clone)
    {
        Name = clone.Name;
        _Case = clone._Case;
        dataStoreOn = clone.dataStoreOn;
        length = clone.length;
        data = clone.data;
    }
    void CopyData(const MemSegment &clone)
    {
        memcpy(data.get(), clone.data.get(), length * sizeof(TField));
    }
    MemSegment DeepCopy()
    {
        MemSegment ret(*this, true);
        return ret;
    }

    // 运算符
    TField &operator[](long i)
    {
        return data[i];
    }
    TField operator[](long i) const
    {
        return data[i];
    }
    MemSegment &operator=(const MemSegment &clone)
    {
        ShallowCopy(clone);
        return *this;
    }
    template <class TOther>
    MemSegment &operator+=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] += right;
        return *this;
    }

    MemSegment &operator+=(const MemSegment &right)
    {
        for (long i = 0; i < length; i++)
            data[i] += right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegment &operator-=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] -= right;
        return *this;
    }

    MemSegment &operator-=(const MemSegment &right)
    {
        for (long i = 0; i < length; i++)
            data[i] -= right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegment &operator*=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] *= right;
        return *this;
    }

    MemSegment &operator*=(const MemSegment &right)
    {
        for (long i = 0; i < length; i++)
            data[i] *= right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegment &operator/=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] /= right;
        return *this;
    }

    MemSegment &operator/=(const MemSegment &right)
    {
        for (long i = 0; i < length; i++)
            data[i] /= right.data[i];
        return *this;
    }

    // 设置值
    void SetValues(TField f)
    {
        for (long i = 0; i < length; i++)
            data[i] = f;
    }

    template <class TOther>
    void SetValues(MemSegment<TCase, TOther> &f)
    {
        for (long i = 0; i < length; i++)
            data[i] = f[i];
    }

    void SetValues(TField *f)
    {
        for (long i = 0; i < length; i++)
            data[i] = f[i];
    }
    void SetValues(TField f, vector<long> &cellIDs)
    {
        size_t n = cellIDs.size();
        for (size_t i = 0; i < n; i++)
        {
            long id = cellIDs[i];
            data[id] = f;
        }
    }
    template <class TOther>
    void AddInplace(TOther f)
    {
        for (long i = 0; i < length; i++)
            data[i] += f;
    }
    template <class TOther>
    void MultiInplace(TOther f)
    {
        for (long i = 0; i < length; i++)
            data[i] *= f;
    }

    template <class TOther>
    void MultiInplaceM(MemSegment<TCase, TOther> &mf)
    {
        for (long i = 0; i < length; i++)
            data[i] *= mf[i];
    }
    template <class TOther>
    void MultiInplaceF(FieldBase<TCase, TOther> *mf)
    {
        if (dataStoreOn == DataStoreOn::onCell)
        {
            for (long i = 0; i < length; i++)
                data[i] *= (*mf->fValues)[i];
        }
        else
        {
            for (long i = 0; i < length; i++)
                data[i] *= (*mf->fFaceValues)[i];
        }
    }
    template <class TOther>
    void DividedInplace(TOther f)
    {
        for (long i = 0; i < length; i++)
            data[i] = f / data[i];
    }
    // 归一化
    void Unit()
    {
        TField maxf = -1e30, minf = 1e-30;
        for (long i = 0; i < length; i++)
        {
            maxf = max(maxf, data[i]);
            minf = min(minf, data[i]);
        }
        auto d = maxf - minf;
        if (abs(d) < 1e-10)
        {
            for (long i = 0; i < length; i++)
            {
                data[i] = (data[i] - minf) / d;
            }
        }
        else
        {
            for (long i = 0; i < length; i++)
            {
                data[i] = (data[i] - minf) / d;
            }
        }
    }
};

// 运算
template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator+(const TOther &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(right);
    tmp += left;
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator+(const MemSegment<TCase, TField> &left, const TOther &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp += right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegment<TCase, TField> operator+(const MemSegment<TCase, TField> &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp += right;
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator-(const TOther &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left - right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator-(const MemSegment<TCase, TField> &left, const TOther &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp -= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegment<TCase, TField> operator-(const MemSegment<TCase, TField> &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp -= right;
    return tmp;
}
template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator*(const TOther &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left * right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator*(const MemSegment<TCase, TField> &left, const TOther &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp *= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegment<TCase, TField> operator*(const MemSegment<TCase, TField> &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp *= right;
    return tmp;
}
template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator/(const TOther &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left / right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegment<TCase, TField> operator/(const MemSegment<TCase, TField> &left, const TOther &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp /= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegment<TCase, TField> operator/(const MemSegment<TCase, TField> &left, const MemSegment<TCase, TField> &right)
{
    MemSegment<TCase, TField> tmp(left);
    tmp /= right;
    return tmp;
}
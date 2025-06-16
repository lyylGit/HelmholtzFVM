#pragma once
#include <iostream>
#include <functional>
#include <cstring>
#include "Mesh3.h"
#include "Types.h"

using namespace std;

#define _NODISCARD [[nodiscard]]

template <class TCase, class TField>
class MemSegmentBC
{
    // 本类主要封装边界面上数据和操作
public:
    // 脚本描述中所使用的符号
    string Name;
    // 关联的Case
    // CaseBase<TCase> *_Case;
    // 数据长度
    long length;
    // 对应边界的起始位置
    long start;
    // 数据
    shared_ptr<TField[]> data;

    MemSegmentBC()
    {
    }
    MemSegmentBC(long _len, long _start)
        : MemSegmentBC("", _len, _start)
    {
    }
    MemSegmentBC(string name, long _len, long _start)
        : Name(name), length(_len), start(_start), data(0)
    {
        data = make_shared<TField[]>(length);
        memset(data.get(), 0, length * sizeof(TField));
    }

    MemSegmentBC(const MemSegmentBC &clone)
    {
        DeepCopy(clone);
    }

    ~MemSegmentBC()
    {
    }
    void SetDimension(long _len, long _start)
    {
        length = _len;
        start = _start;

        data = make_shared<TField[]>(length);
        memset(data.get(), 0, length * sizeof(TField));
    }
    void Zero()
    {
        memset(data.get(), 0, length * sizeof(TField));
    }
    void DeepCopy(const MemSegmentBC &clone)
    {
        Name = clone.Name;
        length = clone.length;
        start =clone. start;
        data = make_shared<TField[]>(length);
        memcpy(data.get(), clone.data.get(), length * sizeof(TField));

        cout << "deepcopy\n";
    }
    void ShallowCopy(const MemSegmentBC &clone)
    {
        Name = clone.Name;
        length = clone.length;
        start =clone. start;
        data = clone.data;
    }
    void CopyData(const MemSegmentBC &clone)
    {
        memcpy(data.get(), clone.data.get(), length * sizeof(TField));
    }
    MemSegmentBC DeepCopy()
    {
        MemSegmentBC ret(*this, true);
        return ret;
    }

    // 运算符
    TField &operator[](long i)
    {
        return data[i - start];
    }
    TField operator[](long i) const
    {
        return data[i - start];
    }
    MemSegmentBC &operator=(const MemSegmentBC &clone)
    {
        ShallowCopy(clone);
        return *this;
    }
    template <class TOther>
    MemSegmentBC &operator+=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] += right;
        return *this;
    }

    MemSegmentBC &operator+=(const MemSegmentBC &right)
    {
        for (long i = 0; i < length; i++)
            data[i] += right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegmentBC &operator-=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] -= right;
        return *this;
    }

    MemSegmentBC &operator-=(const MemSegmentBC &right)
    {
        for (long i = 0; i < length; i++)
            data[i] -= right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegmentBC &operator*=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] *= right;
        return *this;
    }

    MemSegmentBC &operator*=(const MemSegmentBC &right)
    {
        for (long i = 0; i < length; i++)
            data[i] *= right.data[i];
        return *this;
    }
    template <class TOther>
    MemSegmentBC &operator/=(const TOther right)
    {
        for (long i = 0; i < length; i++)
            data[i] /= right;
        return *this;
    }

    MemSegmentBC &operator/=(const MemSegmentBC &right)
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
    void MultiInplace(TOther f)
    {
        for (long i = 0; i < length; i++)
            data[i] *= f;
    }
    template <class TOther>
    void MultiInplaceM(MemSegmentBC<TCase, TOther> mf)
    {
        for (long i = 0; i < length; i++)
            data[i] *= mf[i];
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
_NODISCARD MemSegmentBC<TCase, TField> operator+(const TOther &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(right);
    tmp += left;
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator+(const MemSegmentBC<TCase, TField> &left, const TOther &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp += right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegmentBC<TCase, TField> operator+(const MemSegmentBC<TCase, TField> &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp += right;
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator-(const TOther &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left - right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator-(const MemSegmentBC<TCase, TField> &left, const TOther &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp -= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegmentBC<TCase, TField> operator-(const MemSegmentBC<TCase, TField> &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp -= right;
    return tmp;
}
template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator*(const TOther &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left * right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator*(const MemSegmentBC<TCase, TField> &left, const TOther &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp *= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegmentBC<TCase, TField> operator*(const MemSegmentBC<TCase, TField> &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp *= right;
    return tmp;
}
template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator/(const TOther &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(right);
    for (long i = 0; i < right.length; i++)
        tmp.data[i] = left / right.data[i];
    return tmp;
}

template <class TCase, class TField, class TOther>
_NODISCARD MemSegmentBC<TCase, TField> operator/(const MemSegmentBC<TCase, TField> &left, const TOther &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp /= right;
    return tmp;
}

template <class TCase, class TField>
_NODISCARD MemSegmentBC<TCase, TField> operator/(const MemSegmentBC<TCase, TField> &left, const MemSegmentBC<TCase, TField> &right)
{
    MemSegmentBC<TCase, TField> tmp(left);
    tmp /= right;
    return tmp;
}
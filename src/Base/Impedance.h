#pragma once

#include <memory>
#include "MemSegmentBC.h"

using namespace std;
// 单位面积声阻抗定义
// Z=A0+A1k+A2k^2+A3k^3+...
// k为波数：k=2 * PI * freq / c
template <class TCase, class TField>
class Impedance
{
public:
    int N; // 阶数+1
    // 计算返回的A0,A1,A2,A3....
    shared_ptr<MemSegmentBC<TCase, TField>[]> A;

    // 计算具体的阻抗值，通常与k有关
    function<void(double, Impedance *, void *)> CalcuImpValues;
    // 计算系数，用来计算特征值问题
    function<void(Impedance *)> CalcuImpCoeffs;

public:
    Impedance(int _N) : N(_N)
    {
    }
    Impedance(int _N, int _len, int _start) : N(_N)
    {
        A = make_shared<MemSegmentBC<TCase, TField>[]>(N);
        for (int i = 0; i < _N; i++)
            A[i].SetDimension(_len, _start);
    }
};
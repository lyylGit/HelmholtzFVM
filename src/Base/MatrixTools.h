#pragma once
#include <map>
#include <complex>
#include <string.h>
#include <vector>
#include <fstream>
using namespace std;

// 打印矩阵
template <typename TField>
void _OutputMatrix(map<long, map<long, TField>> &A, TField *B, long n)
{
    ofstream sw("d:/AA.txt", ios::out);
    for (auto it = A.begin(); it != A.end(); it++)
    {
        auto i = it->first;
        auto v = it->second;
        sw << i << "--->\t";
        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            if (abs(itv->second) > 0)
                sw << "【" << itv->first << "】  " << itv->second;
        }

        sw << "<-------------" << B[i] << "\n";
    }
    sw.close();
}
// 统计A的非零个数
template <class TField>
long NNZ(map<long, map<long, TField>> &A)
{
    long nnz = 0;
    for (auto it = A.begin(); it != A.end(); it++)
        nnz += it->second.size();
    return nnz;
}
template <class TField>
void Zero(map<long, map<long, TField>> &A, vector<pair<long, long>> &indexA)
{
    // A.clear();
    // return;

    if (A.empty())
        return;
    if (indexA.empty())
    {
        indexA.reserve(NNZ(A));
        for (auto it = A.begin(); it != A.end(); it++)
        {
            auto i = it->first;
            auto v = it->second;
            for (auto itv = v.begin(); itv != v.end(); itv++)
            {
                A[i][itv->first] = 0;
                indexA.push_back({i, itv->first});
            }
        }
    }
    else
    {
        for (auto ij : indexA)
        {
            A[get<0>(ij)][get<1>(ij)] = 0;
        }
    }
}
template <class TField>
void IndexMatrix(map<long, map<long, TField>> &A, vector<pair<long, long>> &indexA)
{
    if (A.empty())
        return;
    if (indexA.empty())
    {
        indexA.reserve(NNZ(A));
        for (auto it = A.begin(); it != A.end(); it++)
        {
            auto i = it->first;
            auto v = it->second;
            for (auto itv = v.begin(); itv != v.end(); itv++)
            {
                indexA.push_back({i, itv->first});
            }
        }
    }
}
template <class TField>
void Squeeze(map<long, map<long, TField>> &A)
{
    long nnz = 0;

    for (auto it = A.begin(); it != A.end(); it++)
    {
        auto k = it->first;
        auto v = it->second;
        map<long, TField> newA;
        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            int j = itv->first;
            auto p = itv->second;
            if (abs(p) > 1e-30)
            {
                newA[j] = p;
            }
        }
        A[k] = newA;
    }
}
// map 自身是按照key排序
template <class TField>
void FillArrayZ(map<long, map<long, complex<TField>>> &A, double *sa, int *isa, int *jsa, int N)
{
    *(isa) = 1;
    long k = 0;

    for (auto it = A.begin(); it != A.end(); it++)
    {
        int i = it->first;
        auto v = it->second;

        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            int j = itv->first;
            auto p = itv->second;
            *(jsa + k) = j + 1;
            *(sa + 2 * k) = p.real();
            *(sa + 2 * k + 1) = p.imag();
            *(isa + i + 1) = *(isa + i + 1) + 1;
            k++;
        }
    }
    for (int i = 1; i <= N; i++)
    {
        *(isa + i) = *(isa + i) + *(isa + i - 1);
    };
}
template <class TField>
void FillArrayF(map<long, map<long, TField>> &A, double *sa, int *isa, int *jsa, int N, string fn)
{

    *(isa) = 1;
    long k = 0;
    vector<string> L1;
    for (auto it = A.begin(); it != A.end(); it++)
    {

        int i = it->first;
        auto v = it->second;
        string s1 = to_string(i) + "-->\t", s2;
        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            int j = itv->first;
            s1 += "(" + to_string(j) + ")\t";
            auto p = itv->second;
            *(jsa + k) = j + 1;
            *(sa + k) = p.real();
            *(isa + i + 1) = *(isa + i + 1) + 1;

            s1 += to_string(*(sa + k)) + "\t";

            s2 += to_string(*(jsa + k)) + "\t";

            k++;
        }

        L1.push_back(s1 + s2);
    }
    for (int i = 1; i <= N; i++)
    {
        *(isa + i) = *(isa + i) + *(isa + i - 1);
    }
    ofstream sw(fn, ios::out);
    for (int i = 0; i < N; i++)
    {
        sw << L1[i] << "\t" << *(isa + i) << endl;
    }
    sw << *(isa + N) << endl;
    sw.close();
}
template <class TField>
void FillArrayF(map<long, map<long, TField>> &A, double *sa, int *isa, int *jsa, int N)
{
    *(isa) = 1;
    long k = 0;

    for (auto it = A.begin(); it != A.end(); it++)
    {
        int i = it->first;
        auto v = it->second;

        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            int j = itv->first;
            auto p = itv->second;
            *(jsa + k) = j + 1;
            *(sa + k) = p.real();
            *(isa + i + 1) = *(isa + i + 1) + 1;

            k++;
        }
    }
    for (int i = 1; i <= N; i++)
    {
        *(isa + i) = *(isa + i) + *(isa + i - 1);
    };
}

template <class TField>
void WritteArray(map<long, map<long, TField>> &A, int N, string fn)
{
    ofstream sw(fn, ios::out);

    sw << N << "\t" << N << "\t" << NNZ(A) << endl;

    for (auto it = A.begin(); it != A.end(); it++)
    {
        int i = it->first;
        auto v = it->second;

        for (auto itv = v.begin(); itv != v.end(); itv++)
        {
            int j = itv->first;
            auto p = itv->second;

            sw << i + 1 << "\t" << j + 1 << "\t" << p.real() << "\t" << p.imag() << endl;
        }
    }
    sw.close();
}
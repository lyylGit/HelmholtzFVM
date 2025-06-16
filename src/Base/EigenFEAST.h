#pragma once
#include <iostream>
#include <map>
#include <complex>
#include "MatrixTools.h"
#include "Types.h"

using namespace std;

// 右向量 第m mode的第i单元实部
#define getRRe(D, _in, _jm) D[_jm * N * COFACTOR + COFACTOR * _in]
#define getRIm(D, _in, _jm) D[_jm * N * COFACTOR + COFACTOR * _in + 1]
// 左向量第m mode的第i单元实部
#define getLRe(D, _in, _jm) D[M0 * N * COFACTOR + _jm * N * COFACTOR + COFACTOR * _in]
#define getLIm(D, _in, _jm) D[M0 * N * COFACTOR + _jm * N * COFACTOR + COFACTOR * _in + 1]

template <class T>
class Mesh3;
// 计算特征值问题的控制变量
template <class T>
class EigenParam
{
public:
    int N, M0, K;   // 矩阵规模N*N，搜寻特征值数量
    int *fpm;       //[64];
    double *Emid;   //[2];搜寻中心
    double r;       // 搜寻半径
    Mesh3<T> *Mesh; // 输出用

    double *X;       //! eigenvectors
    double *E, *res; //! eigenvalue+residual

    int M, info; // 实际特征值数量，返回信息

    const int COFACTOR = 2; // 复数占2位,实数占1位
    const int LRFACTOR = 2; // factor 2 for L and R, 1 for R only

public:
    EigenParam(int _N, int _M0, int _K, int *_fpm, double *_Emid, double _r, Mesh3<T> *_Mesh)
        : N(_N), M0(_M0), K(_K), fpm(_fpm), Emid(_Emid), r(_r), Mesh(_Mesh)
    {
        for (int i = 0; i < 64; i++)
            fpm[i] = -111;

        E = new double[COFACTOR * M0];                // eigenvalues //factor 2 for complex
        res = new double[M0 * COFACTOR];              // eigenvectors // factor 2 for Left and Right
        X = new double[N * M0 * COFACTOR * LRFACTOR]; // residual //factor 2 for complex // factor 2 for L and R
    }
    ~EigenParam()
    {
        delete[] E;
        delete[] res;
        delete[] X;
    }

    void printInfo()
    {
        /*!!!!!!!!!! REPORT !!!!!!!!!*/
        cout << "FEAST OUTPUT INFO " << info << endl;
        if (info != 0)
            cout << " sparse_zfeast_gcsrgv   -- failed\n";
        if (info == 0)
        {
            double pi2 = 2 * PI;
            cout << "Eigenvalues/Residuals\n";
            for (int i = 0; i <= M - 1; i = i + 1)
            {
                printf("   %d %.15e %.15e %.15e\n", i + 1, *(E + COFACTOR * i) / pi2, *(E + COFACTOR * i + 1) / pi2, *(res + i));
            }
        }
    }
    double *getRightEigenVectors()
    {
        double *D = new double[M * N * COFACTOR];
        for (int j = 0; j < M; j++)
        {
            size_t iMode = j * N * COFACTOR;
            for (long i = 0; i < N; i++)
            {
                double Re = X[iMode + COFACTOR * i], Im = X[iMode + COFACTOR * i + 1];
                D[i * M * COFACTOR + j * COFACTOR] = Re;
                D[i * M * COFACTOR + j * COFACTOR + 1] = Im;
            }
        }
        return D;
    }
    double *getLeftEigenVectors()
    {
        double *D = new double[M * N * COFACTOR];
        for (int j = 0; j < M; j++)
        {
            for (long i = 0; i < N; i++)
            {
                D[i * M * COFACTOR + j * COFACTOR] = getLRe(X, i, j);
                D[i * M * COFACTOR + j * COFACTOR + 1] = getLIm(X, i, j);
            }
        }
        return D;
    }
    void printEigenMatrix(string fn, bool rightVector = true)
    {
        cout << fn << endl;

        ofstream sw(fn, ios::out);
        sw << N << "\t" << M << "\n";

        vector<string> lines;

        if (rightVector)
        {
            for (long i = 0; i < N; i++)
            {
                string s;
                for (int j = 0; j < M; j++)
                {
                    double Re = getRRe(X, i, j),
                           Im = getRIm(X, i, j);

                    s += to_string(Re) + "\t" + to_string(Im) + "\t";
                }
                s += "\n";
                lines.push_back(s);
            }
        }
        else
        {
            for (long i = 0; i < N; i++)
            {
                string s;
                for (int j = 0; j < M; j++)
                {
                    double Re = getLRe(X, i, j),
                           Im = getLIm(X, i, j);

                    s += to_string(Re) + "\t" + to_string(Im) + "\t";
                }
                s += "\n";
                lines.push_back(s);
            }
        }

        for (auto s : lines)
            sw << s;

        sw.close();
    }
    void printEigenMatrixBin(string fn, bool rightVector = true)
    {
        cout << fn << endl;

        ofstream sw(fn, ios::out | ios::binary);
        sw.write(reinterpret_cast<const char *>(&N), sizeof(int));
        sw.write(reinterpret_cast<const char *>(&M), sizeof(int));

        //
        double *D;
        if (rightVector)
            D = getRightEigenVectors();
        else
            D = getLeftEigenVectors();

        sw.write(reinterpret_cast<const char *>(D), M * N * COFACTOR * sizeof(double));
        sw.close();

        delete[] D;
    }
    void printEigenVector(string path, int mode = -1, bool wantLeftVector = false)
    {
        if (mode >= 0)
        {
            printSingleEigenVector(path + "mode_" + to_string(mode), mode, wantLeftVector);
        }
        else
        {
            for (int i = 0; i < M; i++)
            {
                printSingleEigenVector(path + "mode_" + to_string(i), i, wantLeftVector);
            }
        }
    }
    void printSingleEigenVector(string fn, int mode, bool wantLeftVector)
    {
        if (mode >= M)
        {
            cout << "modal:" << mode << " dose not exists\n";
            return;
        }

        double *data = new double[N];

        // R eigenvectors
        for (long i = 0; i < N; i++)
        {
            double Re = getRRe(X, i, mode),
                   Im = getRIm(X, i, mode);

            int m = Re >= 0 ? 1 : -1;

            data[i] = sqrt(Re * Re + Im * Im) * m;
        }
        Mesh->template outputTecPlotCellCenter<double>(fn + "_RA.tec", "pa", data);
        /*  for (long i = 0; i < N; i++)
          {
              data[i] = getRRe(X, i, mode);
          }
          Mesh->template outputTecPlotCellCenter<double>(fn + "_RR.tec", "pa", data);
          for (long i = 0; i < N; i++)
          {
              data[i] = getRIm(X, i, mode);
          }
          Mesh->template outputTecPlotCellCenter<double>(fn + "_RI.tec", "pa", data);
  */
        if (wantLeftVector)
        {
            // L eigenvectors
            for (long i = 0; i < N; i++)
            {
                double Re = getLRe(X, i, mode),
                       Im = getLIm(X, i, mode);
                int m = Re >= 0 ? 1 : -1;
                data[i] = sqrt(Re * Re + Im * Im) * m;
            }
            Mesh->template outputTecPlotCellCenter<double>(fn + "_LA.tec", "pa", data);
         /*   for (long i = 0; i < N; i++)
            {
                data[i] = getLRe(X, i, mode);
            }
            Mesh->template outputTecPlotCellCenter<double>(fn + "_LR.tec", "pa", data);
            for (long i = 0; i < N; i++)
            {
                data[i] = getLIm(X, i, mode);
            }
            Mesh->template outputTecPlotCellCenter<double>(fn + "_LI.tec", "pa", data);
      */  }

         delete[] data;
    }
};
// 计算AX=EBX的特征值,采用FEAST,complex
// 矩阵N*N
template <class TField>
void _GeneralEigenZ(map<long, map<long, TField>> &A, map<long, map<long, TField>> &B, int N);
template <class TField>
void _GeneralEigenF(map<long, map<long, TField>> &A, map<long, map<long, TField>> &B, int N);
template <class TCase, class TField>
void _GeneralEigenZ(map<int, map<long, map<long, TField>>> &A, EigenParam<TCase> &param);

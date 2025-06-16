#pragma once

#define EIGEN_USE_MKL_ALL

#include <Eigen/Sparse>
#include <Eigen/Dense>

#include <fstream>

using namespace std;
using namespace Eigen;

template <typename TField>
void _ConvertMatrix(Eigen::SparseMatrix<TField, RowMajor> &MatrixA, vector<pair<long, long>> &IndexA, map<long, map<long, TField>> &A, long n)
{
    vector<Eigen::Triplet<TField>> tripList;
    tripList.reserve(n * 5);
    for (auto ij : IndexA)
    {
        auto i = get<0>(ij);
        auto j = get<1>(ij);
        tripList.push_back(Eigen::Triplet<TField>(i, j, A[i][j]));
    }

    MatrixA.setFromTriplets(tripList.begin(), tripList.end());
}

template <typename TField>
bool _IterateEquation(map<long, map<long, TField>> &A, vector<pair<long, long>> &IndexA, TField *B, TField *X, long nCells, int nthreads)
{
    if (nthreads > 1)
    {
        Eigen::initParallel();
        Eigen::setNbThreads(nthreads);
    }

    // 系数矩阵
    Eigen::SparseMatrix<TField, RowMajor> MatrixA(nCells, nCells);
    Eigen::Matrix<TField, Eigen::Dynamic, 1> VectorB;

    MatrixA.reserve(VectorXi::Constant(nCells, 6));
    _ConvertMatrix<TField>(MatrixA, IndexA, A, nCells);

    VectorB.resize(nCells);
    for (long i = 0; i < nCells; i++)
        VectorB(i) = B[i];

    // cout << "threads:=" << Eigen::nbThreads() << endl;

    BiCGSTAB<SparseMatrix<TField, RowMajor>, IncompleteLUT<TField>> solver;
    solver.setTolerance(1E-5);
    solver.setMaxIterations(10000);
    solver.compute(MatrixA);
    Eigen::Matrix<TField, Eigen::Dynamic, 1> x = solver.solve(VectorB);
    if (solver.info() != Success)
    {
        cout << "failed! err=" << solver.error() << endl;
        return false;
    }
    else
    {
        cout << "success! err=" << solver.error() << "  n=" << solver.iterations() << endl;
    }
    memcpy(X, x.data(), nCells * sizeof(TField));

    return true;
}
template <typename TField>
bool _IterateEquation(map<long, map<long, TField>> &A, vector<pair<long, long>> &IndexA, TField *B, TField *X, long nCells, TField *C, long M, int nthreads)
{
    if (nthreads > 1)
    {
        Eigen::initParallel();
        Eigen::setNbThreads(nthreads);
    }

    // 系数矩阵
    Eigen::SparseMatrix<TField, RowMajor> MatrixA(nCells, nCells);
    Eigen::Matrix<TField, Eigen::Dynamic, 1> VectorB;

    MatrixA.reserve(VectorXi::Constant(nCells, 6));
    _ConvertMatrix<TField>(MatrixA, IndexA, A, nCells);

    VectorB.resize(nCells);
    for (long i = 0; i < nCells; i++)
        VectorB(i) = B[i];

    // C矩阵
    Map<Eigen::Matrix<TField, Dynamic, Dynamic, RowMajor>> MatrixC(C, nCells, M);

    auto CT = MatrixC.transpose();

    auto tmp = CT * MatrixA * MatrixC;
    Eigen::SparseMatrix<TField, RowMajor> _A = tmp.sparseView();

    auto _B = CT * VectorB;

    BiCGSTAB<SparseMatrix<TField, RowMajor>, IncompleteLUT<TField>> solver;
    solver.setTolerance(1E-7);
    solver.setMaxIterations(10000);
    solver.compute(_A);
    Eigen::Matrix<TField, Eigen::Dynamic, 1> x = solver.solve(_B);
    if (solver.info() != Success)
    {
        cout << "failed! err=" << solver.error() << endl;
        return false;
    }
    else
    {
        cout << "success! err=" << solver.error() << "  n=" << solver.iterations() << endl;
    }
    // 结果返回
    // for (int i = 0; i < M; i++)
    //     x(i) = abs(x(i));
    Eigen::Matrix<TField, Eigen::Dynamic, 1> P = MatrixC * x;
    memcpy(X, P.data(), nCells * sizeof(TField));

    return true;
}
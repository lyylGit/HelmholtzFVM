#include "EigenFEAST.h"
#include <string.h>
#include "feast.h"
#include "feast_sparse.h"
#include <iostream>

// 计算AX=EBX的特征值,采用FEAST,complex
// 矩阵N*N
template <class TField>
void _GeneralEigenZ(map<long, map<long, TField>> &A, map<long, map<long, TField>> &B, int N)
{
    Squeeze(A);
    Squeeze(B);
    // 全部转化成double型
    //  A和B矩阵
    double *sa, *sb;
    int *isa, *jsa, *isb, *jsb;

    int fpm[64];
    double epsout;
    int loop;
    int i, k, err;
    int M0, M, info;
    double Emid[2], r;
    double *X;       //! eigenvectors
    double *E, *res; //! eigenvalue+residual
    long nnz = 0;    // 总的非零单元数

    // A
    nnz = NNZ(A);
    sa = new double[2 * nnz];
    isa = new int[N + 1];
    memset(isa, 0, (N + 1) * sizeof(int));
    jsa = new int[nnz];
    FillArrayZ(A, sa, isa, jsa, N);
    // B
    nnz = NNZ(B);
    sb = new double[2 * nnz];
    isb = new int[N + 1];
    memset(isb, 0, (N + 1) * sizeof(int));
    jsb = new int[nnz];
    FillArrayZ(B, sb, isb, jsb, N);

    /*!!! search interval [Emid,r] including M eigenpairs*/
    Emid[0] = 50;
    Emid[1] = 0.0e0;

    r = 50;
    M0 = 80; // !! M0>=M

    /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
    E = new double[2 * M0];         // eigenvalues //factor 2 for complex
    res = new double[M0 * 2];       // eigenvectors // factor 2 for Left and Right
    X = new double[2 * N * M0 * 2]; // residual //factor 2 for complex // factor 2 for L and R

    /*!!!!!!!!!!!!!!FEAST!!!!!!!!!!!!!*/
    feastinit(fpm);
    fpm[0] = 1; /*change from default value */
    // fpm[7] = 64;
    zfeast_gcsrgv(&N, sa, isa, jsa, sb, isb, jsb, fpm, &epsout, &loop, Emid, &r, &M0, E, X, &M, res, &info);
    cout << "fpm[1]=" << fpm[1] << endl;
    /*!!!!!!!!!! REPORT !!!!!!!!!*/
    cout << "FEAST OUTPUT INFO " << info << endl;
    if (info != 0)
        cout << " sparse_zfeast_gcsrgv   -- failed\n";
    if (info == 0)
    {
        cout << " Csparse_dfeast_gcsrgv   -- success\n";
        cout << "*************************************************\n";
        cout << "************** REPORT ***************************\n";
        cout << "*************************************************\n";
        cout << "Eigenvalues/Residuals\n";
        for (i = 0; i <= M - 1; i = i + 1)
        {
            printf("   %d %.15e %.15e\n", i + 1, *(E + 2 * i), *(E + 2 * i + 1), *(res + i));
        }
    }

    delete[] sa;
    delete[] isa;
    delete[] jsa;

    delete[] sb;
    delete[] isb;
    delete[] jsb;

    delete[] X;
    delete[] E;
    delete[] res;
}

#define _EXPLICIT_DEFINE(CT) template void _GeneralEigenZ(map<long, map<long, CT>> &A, map<long, map<long, CT>> &B, int N)
_EXPLICIT_DEFINE(complex<double>);
_EXPLICIT_DEFINE(complex<float>);
#undef _EXPLICIT_DEFINE

// 矩阵N*N
template <class TField>
void _GeneralEigenF(map<long, map<long, TField>> &A, map<long, map<long, TField>> &B, int N)
{
    Squeeze(A);
    Squeeze(B);
    // 全部转化成double型
    //  A和B矩阵
    double *sa, *sb;
    int *isa, *jsa, *isb, *jsb;

    int fpm[64];
    double epsout;
    int loop;
    int i, k, err;
    int M0, M, info;
    double Emid[2], r;
    double *X;       //! eigenvectors
    double *E, *res; //! eigenvalue+residual
    long nnz = 0;    // 总的非零单元数

    // A
    nnz = NNZ(A);
    sa = new double[nnz];
    isa = new int[N + 1];
    memset(isa, 0, (N + 1) * sizeof(int));
    jsa = new int[nnz];
    FillArrayF(A, sa, isa, jsa, N);
    // FillArrayF(A, sa, isa, jsa, N, "D:\\A.txt");
    // B
    nnz = NNZ(B);
    sb = new double[nnz];
    isb = new int[N + 1];
    memset(isb, 0, (N + 1) * sizeof(int));
    jsb = new int[nnz];
    FillArrayF(B, sb, isb, jsb, N);
    // FillArrayF(B, sb, isb, jsb, N, "D:\\B.txt");

    /*!!! search interval [Emid,r] including M eigenpairs*/
    Emid[0] = 20;
    Emid[1] = 0.0e0;

    r = 50;  // 1e-5;
    M0 = 50; // !! M0>=M

    /*!!!!!!!!!!!!! ALLOCATE VARIABLE */
    E = new double[2 * M0];         // eigenvalues //factor 2 for complex
    res = new double[M0 * 2];       // eigenvectors // factor 2 for Left and Right
    X = new double[2 * N * M0 * 2]; // residual //factor 2 for complex // factor 2 for L and R

    /*!!!!!!!!!!!!!!FEAST!!!!!!!!!!!!!*/
    feastinit(fpm);
    fpm[0] = 1;  /*change from default value */
    fpm[2] = 10; // 收敛误差降低1e-5
    // fpm[7] = 64;
    dfeast_gcsrgv(&N, sa, isa, jsa, sb, isb, jsb, fpm, &epsout, &loop, Emid, &r, &M0, E, X, &M, res, &info);
    // dfeast_gcsrev(&N, sb, isb, jsb, fpm, &epsout, &loop, Emid, &r, &M0, E, X, &M, res, &info);
    // cout << "fpm[1]=" << fpm[1] << endl;
    /*!!!!!!!!!! REPORT !!!!!!!!!*/
    cout << "FEAST OUTPUT INFO " << info << endl;
    if (info != 0)
        cout << " sparse_zfeast_gcsrgv   -- failed\n";
    if (info == 0)
    {
        cout << " Csparse_dfeast_gcsrgv   -- success\n";
        cout << "*************************************************\n";
        cout << "************** REPORT ***************************\n";
        cout << "*************************************************\n";
        cout << "Eigenvalues/Residuals\n";
        for (i = 0; i <= M - 1; i = i + 1)
        {
            printf("   %d %.15e + %.15ei %.15e\n", i + 1, *(E + 2 * i), *(E + 2 * i + 1), *(res + i));
        }
    }

    delete[] sa;
    delete[] isa;
    delete[] jsa;

    delete[] sb;
    delete[] isb;
    delete[] jsb;

    delete[] X;
    delete[] E;
    delete[] res;
}

#define _EXPLICIT_DEFINE(CT) template void _GeneralEigenF(map<long, map<long, CT>> &A, map<long, map<long, CT>> &B, int N)
_EXPLICIT_DEFINE(complex<double>);
_EXPLICIT_DEFINE(complex<float>);
#undef _EXPLICIT_DEFINE

// 计算sigma(λ^i*AiX)=0的特征值,采用FEAST,complex
// 矩阵N*N
template <class TCase, class TField>
void _GeneralEigenZ(map<int, map<long, map<long, TField>>> &A, EigenParam<TCase> &param)
{
    // 全部转化成double型
    //  A和B矩阵
    double *sa;
    int *isa, *jsa;

    double epsout;
    int loop;
    int i, k, err;
   
    int _K = param.K, N = param.N;
    long nnzt = 0; // 最大的非零个数
    for (int i = 0; i < _K; i++)
        nnzt = max(nnzt, NNZ(A[i]));

    sa = new double[_K * 2 * nnzt];
    isa = new int[_K * (N + 1)];
    jsa = new int[_K * nnzt];

    memset(sa, 0, _K * 2 * nnzt * sizeof(double));
    memset(isa, 0, _K * (N + 1) * sizeof(int));
    memset(jsa, 0, _K * nnzt * sizeof(int));

    for (int i = 0; i < _K; i++)
    {
        long nnz = NNZ(A[i]);
        FillArrayZ(A[i], sa + 2 * nnzt * i, isa + (N + 1) * i, jsa + nnzt * i, N);
    }
   
    int _KK = _K - 1;
    char UPLO='F';
    // zfeast_scsrpev(&UPLO,&_KK, &N, sa, isa, jsa, param.fpm, &epsout, &loop, param.Emid, &param.r, &param.M0, param.E, param.X, &param.M, param.res, &param.info);
    // zfeast_hcsrpev(&UPLO,&_KK, &N, sa, isa, jsa, param.fpm, &epsout, &loop, param.Emid, &param.r, &param.M0, param.E, param.X, &param.M, param.res, &param.info);
    zfeast_gcsrpev(&_KK, &N, sa, isa, jsa, param.fpm, &epsout, &loop, param.Emid, &param.r, &param.M0, param.E, param.X, &param.M, param.res, &param.info);

    delete[] sa;
    delete[] isa;
    delete[] jsa;
}

#define _EXPLICIT_DEFINE(CT, ET) template void _GeneralEigenZ(map<int, map<long, map<long, CT>>> &A, EigenParam<ET> &param)
_EXPLICIT_DEFINE(complex<double>, double);
_EXPLICIT_DEFINE(complex<float>, float);
#undef _EXPLICIT_DEFINE
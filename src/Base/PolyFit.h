#pragma once

#include <memory>

using namespace std;
// Adapted from https://github.com/nasa/polyfit/blob/main/PolyFit.cpp
// 拟合多项式 ∑_(i=0)^N(Ai*k^i)
//
template <class T>
class PolyFit
{
public:
    int k;    // 多项式阶数：a0+a1x+a2x^2+...+akx^k
    int n;    // 点数，外部送入
    T *x, *y; // 外部送入

    shared_ptr<T[]> coefbeta; // k+1;Coefficients of the polynomial,Standard error on coefficients
    double SE;                // Standard error
    double R2;

    T **XTWXInv;      // Matrix XTWX Inverse [k+1,k+1]
    double **Weights; // Matrix Weights [n,n]

public:
    PolyFit(int _k, int _n, T *_x, T *_y) : k(_k), n(_n), x(_x), y(_y)
    {
    }
    ~PolyFit()
    {
        Delete2DArray<T>(XTWXInv, k + 1);
        Delete2DArray<double>(Weights, n);
    }
    void Fit()
    {
        // 分配内存
        coefbeta = make_shared<T[]>(k + 1);

        XTWXInv = Make2DArray<T>(k + 1, k + 1);
        Weights = Make2DArray<double>(n, n);

        //
        for (size_t i = 0; i < n; i++)
        {
            Weights[i][i] = 1.;
        }
        // Calculate the coefficients of the fit
        // **************************************************************
        polyFit(x, y, n, k, coefbeta.get(), Weights, XTWXInv);
        double RSS = CalculateRSS(x, y, coefbeta.get(), Weights, n, k + 1);
        double R2 = CalculateR2COD(x, y, coefbeta.get(), Weights, n, k + 1);
        SE = sqrt(RSS / (n - k));

        // Display polynomial
        // **************************************************************
        // DisplayPolynomial(k);

        // Display polynomial coefficients
        // **************************************************************
        // DisplayCoefs(k, coefbeta.get());

        // Display statistics
        // **************************************************************
        // DisplayStatistics(n, k, RSS, R2, SE);
    }

private:
    // Initialize a 2D array
    // **************************************************************
    template <class T2>
    T2 **Make2DArray(const size_t rows, const size_t cols)
    {
        T2 **array;

        array = new T2 *[rows];
        for (size_t i = 0; i < rows; i++)
        {
            array[i] = new T2[cols];
            memset(array[i], 0, cols * sizeof(T2));
        }

        return array;
    }
    template <class T2>
    void Delete2DArray(T2 **&a, int n)
    {
        if (0 != a)
        {
            if (n > 0)
                for (int i = 0; i < n; i++)
                    delete[] (a[i]);
            delete[] a;
        }
        a = 0;
    };

    // Transpose a 2D array
    // **************************************************************
    T **MatTrans(T **array, const size_t rows, const size_t cols)
    {

        T **arrayT = Make2DArray<T>(cols, rows);

        for (size_t i = 0; i < rows; i++)
        {
            for (size_t j = 0; j < cols; j++)
            {
                arrayT[j][i] = array[i][j];
            }
        }

        return arrayT;
    }

    // Perform the multiplication of matrix A[m1,m2] by B[m2,m3]
    // **************************************************************
    template <class T2>
    T **MatMul(const size_t m1, const size_t m2, const size_t m3, T **A, T2 **B)
    {

        T **array = Make2DArray<T>(m1, m3);

        for (size_t i = 0; i < m1; i++)
        {
            for (size_t j = 0; j < m3; j++)
            {
                array[i][j] = 0.;
                for (size_t m = 0; m < m2; m++)
                {
                    array[i][j] += A[i][m] * B[m][j];
                }
            }
        }
        return array;
    }

    // Perform the multiplication of matrix A[m1,m2] by vector v[m2,1]
    // **************************************************************
    void MatVectMul(const size_t m1, const size_t m2, T **A, T *v, T *Av)
    {

        for (size_t i = 0; i < m1; i++)
        {
            Av[i] = 0.;
            for (size_t j = 0; j < m2; j++)
            {
                Av[i] += A[i][j] * v[j];
            }
        }
    }

    // Calculates the determinant of a matrix
    // **************************************************************
    T determinant(T **a, const size_t k)
    {
        T s = 1.0;
        T det = 0.;
        T **b = Make2DArray<T>(k, k);
        size_t m;
        size_t n;

        if (k == 1)
            return (a[0][0]);

        for (size_t c = 0; c < k; c++)
        {

            m = 0;
            n = 0;

            for (size_t i = 0; i < k; i++)
            {

                for (size_t j = 0; j < k; j++)
                {

                    b[i][j] = 0;

                    if (i != 0 && j != c)
                    {

                        b[m][n] = a[i][j];
                        if (n < (k - 2))
                        {
                            n++;
                        }
                        else
                        {
                            n = 0;
                            m++;
                        }
                    }
                }
            }

            det = det + s * (a[0][c] * determinant(b, k - 1));
            s = -1.0 * s;
        }
        Delete2DArray(b, k);
        return (det);
    }

    // Perform the
    // **************************************************************
    void transpose(T **num, T **fac, T **inverse, const size_t r)
    {

        T **b = Make2DArray<T>(r, r);
        T deter;

        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < r; j++)
            {
                b[i][j] = fac[j][i];
            }
        }

        deter = determinant(num, r);

        for (size_t i = 0; i < r; i++)
        {
            for (size_t j = 0; j < r; j++)
            {
                inverse[i][j] = b[i][j] / deter;
            }
        }

        Delete2DArray(b, r);
    }

    // Calculates the cofactors
    // **************************************************************
    void cofactor(T **num, T **inverse, const size_t f)
    {
        T **b = Make2DArray<T>(f, f);
        T **fac = Make2DArray<T>(f, f);

        size_t m;
        size_t n;

        for (size_t q = 0; q < f; q++)
        {

            for (size_t p = 0; p < f; p++)
            {

                m = 0;
                n = 0;

                for (size_t i = 0; i < f; i++)
                {

                    for (size_t j = 0; j < f; j++)
                    {

                        if (i != q && j != p)
                        {

                            b[m][n] = num[i][j];

                            if (n < (f - 2))
                            {
                                n++;
                            }
                            else
                            {
                                n = 0;
                                m++;
                            }
                        }
                    }
                }
                fac[q][p] = pow(-1, q + p) * determinant(b, f - 1);
            }
        }

        transpose(num, fac, inverse, f);

        Delete2DArray(b, f);
        Delete2DArray(fac, f);
    }

    // Display a matrix
    // **************************************************************
    void displayMat(T **A, const size_t n, const size_t m)
    {

        cout << "Matrix " << n << " x " << m << endl;
        for (size_t i = 0; i < n; i++)
        {
            for (size_t j = 0; j < m; j++)
                cout << A[i][j] << "\t";
            cout << endl;
        }
        cout << endl;
    }
    // Perform the fit of data n data points (x,y) with a polynomial of order k
    // **************************************************************
    void polyFit(const T *x, T *y, const size_t n, const size_t k, T *beta, double **Weights, T **XTWXInv)
    {

        // Definition of variables
        // **************************************************************
        T **X = Make2DArray<T>(n, k + 1); // [n,k+1]
        T **XT;                           // [k+1,n]
        T **XTW;                          // [k+1,n]
        T **XTWX;                         // [k+1,k+1]

        T *XTWY = new T[k + 1];
        T *Y = new T[n];

        // Initialize X
        // **************************************************************
        for (size_t i = 0; i < n; i++)
        {
            X[i][0] = 1.0;
            for (size_t j = 1; j < (k + 1); j++)
            {
                X[i][j] = pow(x[i], j);
            }
        }

        // Matrix calculations
        // **************************************************************
        XT = MatTrans(X, n, k + 1);             // Calculate XT
        XTW = MatMul(k + 1, n, n, XT, Weights); // Calculate XT*W
        XTWX = MatMul(k + 1, n, k + 1, XTW, X); // Calculate (XTW)*X

        cofactor(XTWX, XTWXInv, k + 1); // Calculate (XTWX)^-1

        for (size_t m = 0; m < n; m++)
        {
            Y[m] = y[m];
        }
        MatVectMul(k + 1, n, XTW, Y, XTWY);            // Calculate (XTW)*Y
        MatVectMul(k + 1, k + 1, XTWXInv, XTWY, beta); // Calculate beta = (XTWXInv)*XTWY

        // cout << "Matrix X" << endl;
        // displayMat(X, n, k + 1);

        // cout << "Matrix XT" << endl;
        // displayMat(XT, k + 1, n);

        // cout << "Matrix XTW" << endl;
        // displayMat(XTW, k + 1, n);

        // cout << "Matrix XTWXInv" << endl;
        // displayMat(XTWXInv, k + 1, k + 1);

        Delete2DArray(XT, k + 1);
        Delete2DArray(XTW, k + 1);
        Delete2DArray(XTWX, k + 1);
    }
    // Calculate the residual sum of squares (RSS)
    // **************************************************************
    double CalculateRSS(const T *x, const T *y, const T *a, double **Weights, const size_t N, const size_t n)
    {
        T r2 = 0.;
        T ri;
        for (size_t i = 0; i < N; i++)
        {
            ri = y[i];
            ri -= a[0];
            for (size_t j = 1; j < n; j++)
            {
                ri -= a[j] * pow(x[i], j);
            }
            r2 += ri * ri * Weights[i][i];
        }

        return abs(r2);
    }
    // Calculate the total sum of squares (TSS)
    // **************************************************************
    double CalculateTSS(const T *x, const T *y, const T *a, double **Weights, const size_t N, const size_t n)
    {

        T r2 = 0.;
        T ri = 0.;
        T sumwy = 0.;
        T sumweights = 0.;
        size_t begin = 0;
        if (fixed)
        {
            for (size_t i = begin; i < N; i++)
            {
                r2 += y[i] * y[i] * Weights[i][i];
            }
        }
        else
        {

            for (size_t i = begin; i < N; i++)
            {
                sumwy += y[i] * Weights[i][i];
                sumweights += Weights[i][i];
            }

            for (size_t i = begin; i < N; i++)
            {
                ri = y[i] - sumwy / sumweights;
                r2 += ri * ri * Weights[i][i];
            }
        }

        return abs(r2);
    }
    // Calculate coefficient R2 - COD
    // **************************************************************
    double CalculateR2COD(const T *x, const T *y, const T *a, double **Weights, const size_t N, const size_t n)
    {

        double RSS = CalculateRSS(x, y, a, Weights, N, n);
        double TSS = CalculateTSS(x, y, a, Weights, N, n);
        R2 = 1. - RSS / TSS;

        return R2;
    }
    // Display the polynomial
    // **************************************************************
    void DisplayPolynomial(const size_t k)
    {

        cout << "y = ";
        for (size_t i = 0; i < (k + 1); i++)
        {
            cout << "A" << i;
            if (i > 0)
                cout << "X";
            if (i > 1)
                cout << "^" << i;
            if (i < k)
                cout << " + ";
        }
        cout << endl
             << endl;
    }
    // Display the coefficients of the polynomial
    // **************************************************************
    void DisplayCoefs(const size_t k, const T *coefbeta)
    {
        cout << "Polynomial coefficients" << endl;
        cout << "Coeff\tValue" << endl;

        for (size_t i = 0; i < (k + 1); i++)
        {
            cout << "A" << i << "\t";
            cout << coefbeta[i] << "\t";
            cout << endl;
        }
    }
    // Display some statistics values
    // **************************************************************
    void DisplayStatistics(const size_t n, const size_t k, const double RSS, const double R2, const double SE)
    {
        cout << endl;
        cout << "Statistics" << endl;
        cout << "Number of points: " << n << endl;
        cout << "Residual sum of squares: " << RSS << endl;
        cout << "R-square (COD): " << R2 << endl;
        cout << "RMSE: " << SE << endl
             << endl;
    }
};

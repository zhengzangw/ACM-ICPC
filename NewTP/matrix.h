#ifndef MATRIX_H
#define MATRIX_H

#include <bits/stdc++.h>
using namespace std;
#include "math.h"
typedef double MATRIXTYPE;

/**
 * Matrix
 */
struct Matrix {
    MATRIXTYPE m[MAXN][MAXN];
    MATHTYPE l, r;
    Matrix(MATRIXTYPE w) {
        l = r = w;
        memset(m, 0, sizeof(m));
    }
    Matrix() {
        l = r = 0;
        memset(m, 0, sizeof(m));
    }

    void print() {
        for (MATHTYPE i = 0; i < l; ++i, puts("")) {
            printf("%lf", m[i][0]);
            for (MATHTYPE j = 1; j < r; ++j) printf(" %lf", m[i][j]);
        }
    }

    friend Matrix operator*(const Matrix& a, const Matrix& b) {
        Matrix c;
        memset(c.m, 0, sizeof c.m);
        c.l = a.l, c.r = b.r;
        for (MATHTYPE i = 0; i < a.l; ++i)
            for (MATHTYPE j = 0; j < b.r; ++j)
                for (MATHTYPE k = 0; k < a.r; ++k)
                    c.m[i][j] = (c.m[i][j] + a.m[i][k] * a.m[k][j]);
        return c;
    }

    friend Matrix operator*=(Matrix& a, const Matrix& b) {
        a = a * b;
        return a;
    }

    friend Matrix operator+(const Matrix& a, const Matrix& b) {
        assert(a.l == b.l);
        assert(a.r == b.r);
        Matrix c;
        c.l = a.l;
        c.r = a.r;
        for (MATHTYPE i = 0; i < c.l; ++i)
            for (MATHTYPE j = 0; j < c.r; ++j) c.m[i][j] = a.m[i][j] + b.m[i][j];
        return c;
    }

    friend Matrix operator+=(Matrix& a, const Matrix& b) {
        a = a + b;
        return a;
    }

    friend Matrix operator^(const Matrix a, MATHTYPE t) {
        Matrix ans;
        ans.l = ans.r = a.l;
        for (MATHTYPE i = 0; i < a.l; i++) ans.m[i][i] = 1;
        Matrix tmp = a;
        while (t) {
            if (t & 1) ans = tmp * a;
            tmp = tmp * tmp;
            t >>= 1;
        }
        return ans;
    }
};


/**
 * Matrix Inversion
 * @return A^-1
 */
Matrix inv(Matrix a, MATHTYPE n) {
    MATHTYPE m = n << 1;
    for (MATHTYPE i = 0; i < n; i++) a.m[i][n + i] = 1;
    for (MATHTYPE k = 0; k < n; k++) {
        for (MATHTYPE i = k; i < n; i++)
            if (a.m[i][k]) {
                for (MATHTYPE j = 1; j < m; j++) swap(a.m[k][j], a.m[i][j]);
                break;
            }
        if (!a.m[k][k]) {
            puts("No Solution");
            return 0;
        }
        MATHTYPE r = fp(a.m[k][k], mod - 2, mod);
        for (MATHTYPE j = k; j < m; ++j) a.m[k][j] = a.m[k][j] * r;
        for (MATHTYPE i = 0; i < n; i++) {
            r = a.m[i][k];
            if (i != k)
                for (MATHTYPE j = k; j < m; j++)
                    a.m[i][j] = a.m[i][j] - r * a.m[k][j];
        }
    }
    Matrix a_inv(n);
    for (MATHTYPE i = 0; i < n; i++)
        for (MATHTYPE j = 0; j < n; j++) a_inv.m[i][j] = a.m[i][n + j];
    return a_inv;
}


/**
 * The determinant of A
 */
double det(MATHTYPE n, double a[][MAXN]) {
    double maxp;
    MATHTYPE row;
    for (MATHTYPE k = 0; k < n; ++k) {
        maxp = 0;
        for (MATHTYPE i = k; i < n; ++i) {
            if (fabs(a[i][k]) > fabs(maxp)) maxp = a[row = i][k];
        }
        if (fabs(maxp) < eps) return 0;
        if (row != k) {
            for (MATHTYPE j = k; j < n; j++) {
                MATHTYPE t = a[k][j];
                a[k][j] = a[row][j];
                a[row][j] = t;
            }
        }
        for (MATHTYPE j = k + 1; j < n; j++) {
            a[k][j] /= maxp;
            for (MATHTYPE i = k + 1; i < n; i++) a[i][j] -= a[i][k] * a[k][j];
        }
    }

    double ans = 1;
    for (MATHTYPE i = 0; i < n; ++i) ans *= a[i][i];
    return ans;
}


/**
 * Gaussian elimination
 */
bool gauss_cpivot(MATHTYPE n, double a[][MAXN], double b[]) {
    // Solve A*X=b, Ans is in b[] O(n^3)
    MATHTYPE i, j, k, row;
    double maxp, t;
    // For every column
    for (k = 0; k < n; k++) {
        // Find the max element
        for (maxp = 0, i = k; i < n; i++)
            if (fabs(a[i][k]) > fabs(maxp)) maxp = a[row = i][k];
        // If failed, then there is no only solution
        if (fabs(maxp) < eps) return 0;
        // Exchange Row row and Row k
        if (row != k) {
            for (j = k; j < n; j++) {
                t = a[k][j];
                a[k][j] = a[row][j];
                a[row][j] = t;
            }
            t = b[k];
            b[k] = b[row];
            b[row] = t;
        }
        // Update rest elements
        for (j = k + 1; j < n; j++) {
            a[k][j] /= maxp;
            for (i = k + 1; i < n; i++) a[i][j] -= a[i][k] * a[k][j];
        }
        b[k] /= maxp;
        for (i = k + 1; i < n; i++) b[i] -= b[k] * a[i][k];
    }
    // Deal with the 上三角矩阵
    for (i = n - 1; i >= 0; i--)
        for (j = i + 1; j < n; j++) b[i] -= a[i][j] * b[j];
    return 1;
}


#endif
#ifndef POLYNOMIAL_H
#define POLYNOMIAL_H

#include <bits/stdc++.h>
using namespace std;
#include "math.h"


/**
 * Complex
 */
//#if __cplusplus > 199711L
//typedef complex<double> mycomplex;
struct mycomplex {
    double x, y;
    mycomplex(double xx = 0, double yy = 0) : x(xx), y(yy) {}
    friend mycomplex operator+(mycomplex a, mycomplex b) {
        return mycomplex(a.x + b.x, a.y + b.y);
    }
    friend mycomplex operator-(mycomplex a, mycomplex b) {
        return mycomplex(a.x - b.x, a.y - b.y);
    }
    friend mycomplex operator*(mycomplex a, mycomplex b) {
        return mycomplex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
};


/**
 * NTT
 * @timecomplex O(nlgn)
 * @const 998244353 = 119*2^23 + 1; g = 3;gni2 = 332748118;
 * @const 1004535809 = 479*2^21 + 1; g = 3; gni3 = 334845270;
 * @const 2281701377 = 17*2^27 + 1; g = 3;
 * @const 469762049 = 7*2^26 + 1; g = 3; gni1 = 156587350;
 */

// For MTT
/*
MATHTYPE Pnum = 1;
MATHTYPE Proot_cur,Proot_invcur,Pcur;
void set_prime(){
    Proot_cur=Proot[Pnum];
    Proot_invcur=Proot_inv[Pnum];
    Pcur=P[Pnum];
}
*/
MATHTYPE r[MAXN * 2];
void calrev(MATHTYPE n, MATHTYPE len) {
    for (MATHTYPE i = 0; i < n; ++i)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (len - 1));
}
const MATHTYPE P[] = {469762049, 998244353, 1004535809};
const MATHTYPE Proot[] = {3, 3, 3};
const MATHTYPE Proot_inv[] = {156587350, 332748118, 334845270};
const MATHTYPE Pnum = 1;
const MATHTYPE Proot_cur = Proot[Pnum], Proot_invcur = Proot_inv[Pnum],
         Pcur = P[Pnum];

void NTT(MATHTYPE limit, MATHTYPE* A, MATHTYPE type) {
    for (MATHTYPE i = 0; i < limit; ++i)
        if (i < r[i]) std::swap(A[i], A[r[i]]);
    for (MATHTYPE mid = 1; mid < limit; mid <<= 1) {
        MATHTYPE wn = fp(type == 1 ? Proot_cur : Proot_invcur,
                   (Pcur - 1) / (mid << 1), Pcur);
        for (MATHTYPE j = 0; j < limit; j += (mid << 1)) {
            MATHTYPE w = 1;
            for (MATHTYPE k = 0; k < mid; k++, w = (w * wn) % Pcur) {
                MATHTYPE x = A[j + k], y = w * A[j + k + mid] % Pcur;
                A[j + k] = (x + y) % Pcur;
                A[j + k + mid] = (x - y + Pcur) % Pcur;
            }
        }
    }
    if (type == -1) {
        MATHTYPE inv = fp(limit, Pcur - 2, Pcur);
        for (MATHTYPE i = 0; i < limit; ++i) A[i] = (A[i] * inv) % Pcur;
    }
}


/**
 * polynomial multiplation
 */
MATHTYPE sa[MAXN], sb[MAXN];
void poly_mult(MATHTYPE n, MATHTYPE m, MATHTYPE a[], MATHTYPE b[], MATHTYPE ans[]) {
    memset(sa, 0, sizeof(sa));
    memset(sb, 0, sizeof(sb));
    // set_prime();
    MATHTYPE limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    for (MATHTYPE i = 0; i < limit; ++i) {
        sa[i] = a[i];
        sb[i] = b[i];
    }
    calrev(limit, len);
    NTT(limit, sa, 1);
    NTT(limit, sb, 1);
    for (MATHTYPE i = 0; i < limit; i++) ans[i] = (sa[i] * sb[i]) % Pcur;
    NTT(limit, ans, -1);
}


/**
 * polynomial inversion
 * @return G(x)*F(x) = 1 mod x^n
 */
MATHTYPE t1[MAXN], t2[MAXN];
void poly_inv(MATHTYPE* a, MATHTYPE n, MATHTYPE mod, MATHTYPE* ans) {
    memset(t1, 0, sizeof(t1));
    memset(t2, 0, sizeof(t2));
    // set_prime();
    ans[0] = fp(a[0], mod - 2, mod);
    for (MATHTYPE bas = 2, len = 2, limit = 4; bas < (n << 1);
         bas <<= 1, ++len, limit <<= 1) {
        calrev(limit, len);
        for (MATHTYPE i = 0; i < bas; ++i) {
            t1[i] = a[i];
            t2[i] = ans[i];
        }
        NTT(limit, t1, 1);
        NTT(limit, t2, 1);
        for (MATHTYPE i = 0; i < limit; ++i)
            ans[i] = md(t2[i] * md(2 - md(t1[i] * t2[i])));
        NTT(limit, ans, -1);
        for (MATHTYPE i = 0; i < bas; ++i) ans[i] %= mod;
        for (MATHTYPE i = bas; i < limit; ++i) ans[i] = 0;
    }
}


/**
 * polynomial division
 * @return F(x) = Q(x)G(x) + R(x) n = m*(n-m) + (m-1)
 * @timecomplex O(nlgn)
 */
MATHTYPE G_R[MAXN], G_R_inv[MAXN], F_R[MAXN], Q_R[MAXN], G_Q[MAXN];
void poly_div(MATHTYPE n, MATHTYPE m, MATHTYPE mod, MATHTYPE F[], MATHTYPE G[], MATHTYPE Q[], MATHTYPE R[]) {
    for (MATHTYPE i = 0; i <= m; ++i) G_R[i] = G[m - i];
    for (MATHTYPE i = n - m + 1; i <= n; ++i) G_R[i] = 0;
    poly_inv(G_R, n - m + 1, mod, G_R_inv);
    for (MATHTYPE i = 0; i <= n - m; ++i) F_R[i] = F[n - i];
    poly_mult(n - m, n - m, F_R, G_R_inv, Q_R);
    for (MATHTYPE i = 0; i <= n - m; ++i) Q[i] = Q_R[n - m - i];
    poly_mult(m, n - m, G, Q, G_Q);
    for (MATHTYPE i = 0; i <= n; ++i) R[i] = md(F[i] - G_Q[i]);
}


/**
 * MTT
 * NTT for any mod rather than prime
 */
const MATHTYPE P1 = P[0], P2 = P[1], P3 = P[2];
const MATHTYPE M = P1 * P2;
MATHTYPE a1[MAXN], a2[MAXN], a3[MAXN], b1[MAXN], b2[MAXN], b3[MAXN];
void MTT(MATHTYPE n, MATHTYPE m, MATHTYPE mod, MATHTYPE a[], MATHTYPE b[], MATHTYPE ans[]) {
    for (MATHTYPE i = 0; i <= n; ++i) a1[i] = a2[i] = a3[i] = a[i];
    for (MATHTYPE i = 0; i <= m; ++i) b1[i] = b2[i] = b3[i] = b[i];
    MATHTYPE limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    calrev(limit, len);
    // Pnum = 0;
    poly_mult(n, m, a1, b1, a1);
    // Pnum = 1;
    poly_mult(n, m, a2, b2, a2);
    // Pnum = 2;
    poly_mult(n, m, a3, b3, a3);
    for (MATHTYPE i = 0; i <= n + m; ++i) {
        MATHTYPE tmp = (multiMod((P2 * a1[i]), inv(P2, P1), M) +
                  multiMod((P1 * a2[i]) % M, inv(P1, P2), M)) %
                 M;
        MATHTYPE k = ((a3[i] - tmp) % P3 + P3) % P3 * inv(M, P3) % P3;
        ans[i] = ((k % mod) * (M % mod) % mod + tmp % mod) % mod;
    }
}


/**
 * FFT
 * @return C_k = \sun_{i+j=k}(A_i+A_j)
 * @timecomplex O(nlgn)
 * -Way1 iteration
 * -(Deprecated) divide and conquer, A(w)=A'(w^2)+wA(w^2) A(w+n/2) = A'(w^2)-wA(w^2)
 */
void FFT(MATHTYPE limit, mycomplex* a, MATHTYPE type) {
    mycomplex x, y;
    for (MATHTYPE i = 0; i < limit; ++i)
        if (i < r[i]) std::swap(a[i], a[r[i]]);
    for (MATHTYPE mid = 1; mid < limit; mid <<= 1) {
        mycomplex Wn(cos(pi / mid), type * sin(pi / mid)), x, y;
        for (MATHTYPE j = 0; j < limit; j += (mid << 1)) {
            mycomplex w(1, 0);
            for (MATHTYPE k = 0; k < mid; ++k, w = w * Wn) {
                x = a[j + k];
                y = w * a[j + mid + k];
                a[j + k] = x + y;
                a[j + mid + k] = x - y;
            }
        }
    }
}


/**
 * Polynomial Multiplication
 */
mycomplex ta[MAXN], tb[MAXN], tc[MAXN];
void poly_mult_double(MATHTYPE n, MATHTYPE m, MATHTYPE F[], MATHTYPE G[], MATHTYPE ans[]) {
    MATHTYPE limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    for (MATHTYPE i = 0; i <= limit; ++i) {
        ta[i] = F[i];
        tb[i] = G[i];
        tc[i] = 0;
    }
    //二进制反转
    calrev(limit, len);
    FFT(limit, ta, 1);
    FFT(limit, tb, 1);
    for (MATHTYPE i = 0; i < limit; ++i) tc[i] = ta[i] * tb[i];
    FFT(limit, tc, -1);
    for (MATHTYPE i = 0; i <= n + m; ++i) ans[i] = (MATHTYPE)(tc[i].x / limit + 0.5);
}


/**
 * Polynomial Module
 */
MATHTYPE TG[MAXN], TR[MAXN];
void poly_mult_mod(MATHTYPE m, MATHTYPE* A, MATHTYPE* B, MATHTYPE* MOD, MATHTYPE opt) {
    opt = min(opt, m << 1);
    memset(TG, 0, sizeof(TG));
    memset(TR, 0, sizeof(TR));
    poly_mult(m, m, A, B, A);
    if (opt > m) {
        poly_div(opt, m, mod, A, MOD, TG, TR);
        for (MATHTYPE i = m; i <= m << 1; ++i) A[i] = 0;
        for (MATHTYPE i = 0; i < m; ++i) A[i] = TR[i];
    }
}


/**
 * Polynomial Fast power
 * @return F(x)^n % G(x)
 */
void poly_fast_pow(MATHTYPE* A, MATHTYPE n, MATHTYPE* ans, MATHTYPE* MOD, MATHTYPE m) {
    for (MATHTYPE i = 1; n; n >>= 1, poly_mult_mod(m, A, A, MOD, 2 << i), ++i) {
        if (n & 1) poly_mult_mod(m, ans, A, MOD, m << 1);
        if ((1 << (i - 1)) >= (m << 1)) --i;
    }
}


/**
 * Solve Linear Recursive Equation
 * @timecomplex O(mlgmlgn)
 */
MATHTYPE T[MAXN], C[MAXN], G[MAXN];
MATHTYPE Solve_Linear_Recur(MATHTYPE n, MATHTYPE m, MATHTYPE* f, MATHTYPE* A) {
    for (MATHTYPE i = 1; i <= m; ++i) G[m - i] = (mod - f[i]) % mod;
    G[m] = 1;
    T[1] = 1;
    C[0] = 1;
    poly_fast_pow(T, n, C, G, m);
    MATHTYPE ans = 0;
    for (MATHTYPE i = 0; i <= m - 1; ++i) ans = (ans + md(A[i] * C[i])) % mod;
    return ans;
}


/**
 * Fast Walsh-Hadamard Transform
 * @timecomplex O(nlgn)
 */
inline MATHTYPE div2(MATHTYPE x) { return x & 1 ? (mod + x) / 2 : x / 2; }
void FWTxor(MATHTYPE limit, MATHTYPE* a, MATHTYPE type) {
    MATHTYPE t;
    for (MATHTYPE mid = 1; mid < limit; mid <<= 1)
        for (MATHTYPE j = 0; j < limit; j += (mid << 1))
            for (MATHTYPE k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] = type == 1 ? md(t + a[j + k + mid])
                                     : div2(md(t + a[j + k + mid]));
                a[j + k + mid] = type == 1 ? md(t - a[j + k + mid])
                                           : div2(md(t - a[j + k + mid]));
            }
}
void FWTand(MATHTYPE limit, MATHTYPE* a, MATHTYPE type) {
    MATHTYPE t;
    for (MATHTYPE mid = 1; mid < limit; mid <<= 1)
        for (MATHTYPE j = 0; j < limit; j += (mid << 1))
            for (MATHTYPE k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] =
                    type == 1 ? md(t + a[j + k + mid]) : md(t - a[j + k + mid]);
                a[j + k + mid] = md(a[j + k + mid]);
            }
}
void FWTor(MATHTYPE limit, MATHTYPE* a, MATHTYPE type) {
    MATHTYPE t;
    for (MATHTYPE mid = 1; mid < limit; mid <<= 1)
        for (MATHTYPE j = 0; j < limit; j += (mid << 1))
            for (MATHTYPE k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] = md(t);
                a[j + k + mid] =
                    type == 1 ? md(t + a[j + k + mid]) : md(a[j + k + mid] - t);
            }
}


/**
 * Fast Mobius Transform
 * O(nlgn)
 */
void FMTor(MATHTYPE n, MATHTYPE* a, MATHTYPE type) {
    // n = log(limit)
    for (MATHTYPE i = 0; i < n; ++i)
        for (MATHTYPE j = 0; j < 1 << n; j++)
            if (j >> i & 1)
                a[j] = md(a[j] + a[j ^ (1 << i)] * (type == 1 ? 1 : -1));
}
void FMTand(MATHTYPE n, MATHTYPE* a, MATHTYPE type) {
    for (MATHTYPE i = 0; i < n; ++i)
        for (MATHTYPE j = (1 << n) - 1; ~j; j--)
            if (~j >> i & 1)
                a[j] = md(a[j] + a[j | 1 << i] * (type == 1 ? 1 : -1));
}

#endif
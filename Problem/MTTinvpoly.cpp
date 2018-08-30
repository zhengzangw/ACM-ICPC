#include <bits/stdc++.h>
#define rep(i, a, b) for (i = (a); i <= (b); ++i)
#define repi(i, a, b) for (i = (a); i >= (b); --i)
#define eps 1e-10
#define randint(x, y) rand() % (y) + (x)
typedef long long LL;
using namespace std;
const LL mod = 1e9 + 7;
#define md(x) (((x) % mod + mod) % mod)
LL i, j, k, n, m, ans;
const LL MAXN = 5e5 + 10;

inline void read(LL &n) //快读
{
    LL x(0), f(1);
    char c;
    for (c = getchar(); !isdigit(c); c = getchar())
        if (c == '-')
            f = -1;
    for (; isdigit(c); c = getchar())
        x = (x << 3) + (x << 1) + (c ^ 48);
    n = x * f;
}

inline void write(LL u, char ed = '\n')
{
    if (!u)
    {
        putchar(48);
        putchar(ed);
        return;
    }
    static int sta[43], tp;
    for (tp = 0; u; u /= 10)
        sta[++tp] = u % 10;
    for (; tp; putchar(sta[tp--] ^ 48))
        ;
    putchar(ed);
}

//快速幂
LL FM(LL a, LL b, LL mod)
{
    LL c = 1;
    for (; b; b >>= 1, a = 1ll * a * a % mod)
        if (b & 1)
            c = 1ll * c * a % mod;
    return c % mod;
}

//快乘
LL multiMod(LL a, LL b, LL mod)
{
    b %= mod;
    a %= mod;
    LL res = 0;
    while (b)
    {
        if (b & 1)
        {
            res += a;
            if (res >= mod)
                res -= mod;
        }
        b >>= 1;
        a = a << 1;
        if (a >= mod)
            a -= mod;
    }
    return res;
}

//扩展欧几里得
LL ext_gcd(LL a, LL b, LL &x, LL &y)
{
    if (b == 0)
    {
        x = 1;
        y = 0;
        return a;
    }
    else
    {
        LL r = ext_gcd(b, a % b, y, x);
        y -= x * (a / b);
        return r;
    }
}

//逆元
//扩展欧几里得
LL inv(LL a, LL mod)
{ // only for (x,mod)=1
    a %= mod;
    LL x, y;
    ext_gcd(mod, a, x, y);
    return ((y % mod) + mod) % mod;
}

LL r[MAXN * 2];
void calrev(LL n, LL len)
{
    for (LL i = 0; i < n; ++i)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (len - 1));
}

void NTT(LL limit, LL *A, LL type, LL P = 998244353, LL G = 3, LL Gni = 332748118)
{
    for (LL i = 0; i < limit; ++i)
        if (i < r[i])
            std::swap(A[i], A[r[i]]);
    for (LL mid = 1; mid < limit; mid <<= 1)
    {
        LL wn = FM(type == 1 ? G : Gni, (P - 1) / (mid << 1), P);
        for (LL j = 0; j < limit; j += (mid << 1))
        {
            LL w = 1;
            for (LL k = 0; k < mid; k++, w = (w * wn) % P)
            {
                LL x = A[j + k], y = w * A[j + k + mid] % P;
                A[j + k] = (x + y) % P;
                A[j + k + mid] = (x - y + P) % P;
            }
        }
    }
}
void poly_mult(LL n, LL m, LL a[], LL b[], LL P = 998244353, LL G = 3, LL Gni = 332748118)
{
    LL limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    calrev(limit, len);
    NTT(limit, a, 1, P, G, Gni);
    NTT(limit, b, 1, P, G, Gni);
    for (LL i = 0; i < limit; i++)
        a[i] = multiMod(b[i], a[i], P);
    NTT(limit, a, -1, P, G, Gni);
    LL inv = FM(limit, P - 2, P);
    for (LL i = 0; i <= n + m; ++i)
        a[i] = multiMod(a[i], inv, P);
}
//任意模数NTT
LL a1[MAXN], a2[MAXN], a3[MAXN], b1[MAXN], b2[MAXN], b3[MAXN];
void MTT(LL n, LL m, LL F[], LL G[], LL A[], LL mod)
{
    const LL P1 = 469762049, g1 = 3, gni1 = 156587350;
    const LL P2 = 998244353, g2 = 3, gni2 = 332748118;
    const LL P3 = 1004535809, g3 = 3, gni3 = 334845270;
    LL limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    calrev(limit, len);
    for (LL i = n; i <= limit; ++i)
        a1[i] = a2[i] = a3[i] = b1[i] = b2[i] = b3[i] = 0;
    for (LL i = 0; i <= n; ++i)
        a1[i] = a2[i] = a3[i] = F[i];
    for (LL i = 0; i <= m; ++i)
        b1[i] = b2[i] = b3[i] = G[i];
    poly_mult(n, m, a1, b1, P1, g1, gni1);
    poly_mult(n, m, a2, b2, P2, g2, gni2);
    poly_mult(n, m, a3, b3, P3, g3, gni3);
    for (LL i = 0; i <= n + m; ++i)
    {
        LL M = P1 * P2;
        LL ans = (multiMod((P2 * a1[i]), inv(P2, P1), M) + multiMod((P1 * a2[i]) % M, inv(P1, P2), M)) % M;
        LL k = ((a3[i] - ans) % P3 + P3) % P3 * inv(M, P3) % P3;
        A[i] = ((k % mod) * (M % mod) % mod + ans % mod) % mod;
    }
}

//多项式求逆
LL t1[MAXN], t2[MAXN], t3[MAXN], t4[MAXN], deg;

LL a[MAXN], b[MAXN], orz[MAXN];
int main()
{
    read(n);
    rep(i, 0, n - 1) read(a[i]);
    for (deg = 1; deg < n; deg <<= 1)
        ;
    b[0] = FM(a[0], mod - 2, mod);

    for (LL bas = 2, len = 2; bas <= deg; bas <<= 1, ++len)
    {
        for (LL i = 0; i < bas; ++i)
        {
            t1[i] = a[i];
            t2[i] = b[i];
        }
        //b[i] * (((2-a[i]*b[i])%P+P)%P) % P;
        MTT(bas, bas, t1, t2, t3, mod);
        MTT(bas, bas, t2, t3, t4, mod);
        rep(i, 0, bas - 1) b[i] = md((b[i] << 1) - t4[i]);
        //rep(i,0,bas-1) b[i] = md(t3[i]);
    }
    rep(i, 0, n - 1) write(b[i], ' ');
    return 0;
}

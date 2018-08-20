#include <bits/stdc++.h>
using namespace std;
const int MAXN = 3e5 + 10;
typedef long long LL;
LL limit, len;
LL r[2 * MAXN];
LL n;
long long g[MAXN], f[MAXN];
const LL P = 998244353, G = 3, Gni = 332748118;
#define md(x) (((x) % P + P) % P)
LL FM(LL a, LL b, LL mod) //快速幂
{
    LL c = 1;
    for (; b; b >>= 1, a = 1ll * a * a % mod)
        if (b & 1)
            c = 1ll * c * a % mod;
    return c;
}
void NTT(LL limit, LL *A, LL type)
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
void polymulti(LL F[], LL G[])
{
    for (LL i = 0; i < limit; ++i)
        r[i] = (r[i >> 1] >> 1 | ((i & 1) << (len - 1)));
    NTT(limit, F, 1);
    NTT(limit, G, 1);
    for (LL i = 0; i < limit; i++)
        F[i] = (F[i] * G[i]) % P;
    NTT(limit, F, -1);
}

LL a[MAXN], b[MAXN];
void conquer(LL l, LL r)
{
    if (l == r)
        return;
    LL mid = (l + r) >> 1;
    conquer(l, mid);
    for (limit = 1, len = 0; limit <= (r - l + 1) << 1; limit <<= 1, len++)
        ;
    for (LL i = 0; i <= limit; ++i)
        a[i] = b[i] = 0;
    for (LL i = l; i <= mid; ++i)
        a[i - l] = f[i];
    for (LL i = 0; i <= r - l; ++i)
        b[i] = g[i];
    polymulti(a, b);
    LL inv = FM(limit, P - 2, P);
    for (LL i = mid + 1; i <= r; ++i)
    {
        f[i] = md(f[i] + md(a[i - l] * inv));
    }
    conquer(mid + 1, r);
}
int main()
{
    scanf("%lld", &n);
    f[0] = 1;
    g[0] = 0;
    for (LL i = 1; i < n; ++i)
        scanf("%lld", &g[i]);
    conquer(0, n - 1);
    for (LL i = 0; i < n; ++i)
    {
        printf("%lld ", f[i]);
    }
    return 0;
}

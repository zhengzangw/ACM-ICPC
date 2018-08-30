#include <bits/stdc++.h>
#define rep(i, a, b) for (i = (a); i <= (b); ++i)
#define repi(i, a, b) for (i = (a); i >= (b); --i)
#define eps 1e-10
#define randint(x, y) rand() % (y) + (x)
typedef long long LL;
using namespace std;
LL mod;
#define md(x) (((x) % mod + mod) % mod)
LL i,j,k,n,m,ans,p;

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

LL ni(LL a, LL mod)
{ // only for (x,mod)=1
    a %= mod;
    LL x, y;
    ext_gcd(mod, a, x, y);
    return ((y % mod) + mod) % mod;
}

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

inline void write(LL u,char ed = '\n')
{
    if (!u){putchar(48);putchar(ed);return;}
    static int sta[43],tp;
    for (tp=0;u;u/=10) sta[++tp]=u%10;
    for (;tp;putchar(sta[tp--]^48));
    putchar(ed);
}

const LL MAXN = 10e5+10;
LL a1[MAXN],b1[MAXN],a2[MAXN],b2[MAXN],a3[MAXN],b3[MAXN];
LL r[MAXN * 2];
LL FM(LL a,LL b,LL mod){
    LL c = 1;
    for (;b;b>>=1,a=1ll*a*a%mod)
        if (b&1)
            c = 1ll*a*c%mod;
    return c;
}

void NTT(LL limit, LL *A, LL type,LL P,LL G,LL Gni)
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

LL limit,len;
void polymulti(LL n,LL m,LL F[],LL G[],LL P,LL g,LL Gni)
{
    NTT(limit, F, 1, P, g, Gni);
    NTT(limit, G, 1, P, g, Gni);
    for (LL i = 0; i < limit; i++)
        F[i] = (F[i] * G[i]) % P;
    NTT(limit,F, -1, P, g, Gni);
    LL inv = FM(limit, P - 2, P);
    for (LL i = 0; i <= n + m; ++i)
        F[i] = (F[i] * inv) % P;
}
const LL P1 = 469762049, g1 = 3, gni1 = ni(g1,P1);
const LL P2 = 998244353 , g2 = 3, gni2 = ni(g2,P2);
const LL P3 = 1004535809, g3 = 3, gni3 = ni(g3,P3);
LL multiMod(LL a, LL b, LL mod) //快速乘
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

LL MTT(LL a1,LL a2, LL a3){
    LL M = P1*P2;
    LL ans = (multiMod((P2*a1),ni(P2,P1),M)+multiMod((P1*a2)%M,ni(P1,P2),M))%M;
    LL k = ((a3-ans)%P3+P3)%P3 * ni(M,P3) % P3;
    return ((k%mod)*(M%mod) % mod + ans % mod)%mod;
}
int main(){
    read(n); read(m); read(mod);
    rep(i,0,n) read(a1[i]),a2[i]=a1[i],a3[i]=a1[i];
    rep(i,0,m) read(b1[i]),b2[i]=b1[i],b3[i]=b1[i];
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    for (LL i = 0; i < limit; ++i)
        r[i] = (r[i >> 1] >> 1 | ((i & 1) << (len - 1)));
    polymulti(n,m,a1,b1,P1,g1,gni1);
    polymulti(n,m,a2,b2,P2,g2,gni2);
    polymulti(n,m,a3,b3,P3,g3,gni3);
    rep(i,0,n+m) write(MTT(a1[i],a2[i],a3[i]),' ');
    return 0;
}

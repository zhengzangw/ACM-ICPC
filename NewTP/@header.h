#ifndef header_H
#define header_H

#include <bits/stdc++.h>
#define rep(i, a, b) for (i = (a); i <= (b); ++i)
#define repi(i, a, b) for (i = (a); i >= (b); --i)
#define eps 1e-10
#define randint(x, y) rand() % (y) + (x)
typedef long long LL;
using namespace std;
const int mod = 998244353;
#define md(x) (((x) % mod + mod) % mod)
int i,j,k,n,m,ans;

void file(char s[])
{
#ifndef ONLINE_JUDGE
    freopen(strcat(s,".in"), "r", stdin);
    freopen(strcat(s,"out"), "w", stdout);
#endif
}
inline void read(int &n) //快读
{
    int x(0), f(1);
    char c;
    for (c = getchar(); !isdigit(c); c = getchar())
        if (c == '-')
            f = -1;
    for (; isdigit(c); c = getchar())
        x = (x << 3) + (x << 1) + (c ^ 48);
    n = x * f;
}

inline void write(int u,char ed = '\n')
{
    if (!u){putchar(48);putchar(ed);return;}
    static int sta[43],tp;
    for (tp=0;u;u/=10) sta[++tp]=u%10;
    for (;tp;putchar(sta[tp--]^48));
    putchar(ed);
}

#endif
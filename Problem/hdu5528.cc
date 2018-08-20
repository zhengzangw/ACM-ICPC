#include <cstdio>
#include <cmath>
#include <iostream>
const int MAXP = (int)10e5;
int prime[MAXP + 10], p[MAXP + 10];

void pre_ptable()
{
    p[1] = 1;
    int cnt = 0;
    for (int i = 2; i <= MAXP; ++i)
    {
        if (!p[i])
            prime[cnt++] = i;
        for (int j = 0;; ++j)
        {
            int x = i * prime[j];
            if (x > MAXP)
                break;
            p[x] = 1;
            if (i % prime[j] == 0)
                break;
        }
    }
}

unsigned long long split(int n, int a[], int &len)
{
    int sqrtn = (int)sqrt(n);
    len = 0;
    a[0] = 1;
    unsigned long long ans = 1, tmp = 1;
    for (int i = 0; prime[i] <= sqrtn; ++i)
    {
        tmp = 0;
        while (n % prime[i] == 0)
        {
            tmp++;
            a[++len] = prime[i];
            n /= prime[i];
        }
        ans = ans * (tmp + 1);
    }
    if (n != 1)
    {
        a[++len] = n;
        ans = 2 * ans;
    }
    return ans;
}

unsigned long long P(int fac[], int &len)
{
    unsigned long long rnt = 1, tmp = 0, cnt = 1;
    for (int i = 1; i <= len; ++i)
    {
        if (fac[i] != fac[i - 1])
        {
            rnt = rnt * (tmp + 1);
            cnt = fac[i];
            tmp = 0;
        }
        tmp = tmp + cnt * cnt;
        cnt = cnt * fac[i];
    }
    rnt = rnt * (tmp + 1);
    return rnt;
}

int main()
{
    //std::cout << mod;
    int T, n, fac[100000], len;
    unsigned long long ans, Q, Pt;
    pre_ptable();
    scanf("%d", &T);
    while (T--)
    {
        scanf("%d", &n);
        Q = n * split(n, fac, len);
        Pt = P(fac, len);
        ans = P(fac, len) - Q;
        printf("%llu\n", ans);
    }
    return 0;
}

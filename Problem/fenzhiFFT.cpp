#include <cstdio>
#include <algorithm>
#include <cmath>
const int mod = 998244353;
const int MAXN = 10e5+10;
int n,limit;
#define pi acos(-1)
#define md(x) (((x)%mod)+mod)%mod
//复数
struct complex
{
    double x, y;
    complex(double xx = 0, double yy = 0) : x(xx), y(yy) {}
};
complex operator+(complex a, complex b) { return complex(a.x + b.x, a.y + b.y); }
complex operator-(complex a, complex b) { return complex(a.x - b.x, a.y - b.y); }
complex operator*(complex a, complex b) { return complex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
complex g[MAXN], f[MAXN], a[MAXN], b[MAXN];
//FFT O(nlogn) 多项式乘法，两个数组卷积
//分治法慢 A(w)=A'(w^2)+wA(w^2) A(w+n/2) = A'(w^2)-wA(w^2)
//迭代实现 10e6
//卷积 ck = \sum_{i+j=k}(ai+aj)
int len(0), r[MAXN * 2];
void FFT(int limit, complex *a, int type)
{
    complex x, y;
    for (int i = 0; i < limit; ++i)
        if (i < r[i])
            std::swap(a[i], a[r[i]]);
    for (int mid = 1; mid < limit; mid <<= 1)
    {
        complex Wn(cos(pi / mid), type * sin(pi / mid)), x, y;
        for (int j = 0; j < limit; j += (mid << 1))
        {
            complex w(1, 0);
            for (int k = 0; k < mid; ++k, w = w * Wn)
            {
                x = a[j + k];
                y = w * a[j + mid + k];
                a[j + k] = x + y;
                a[j + mid + k] = x - y;
            }
        }
    }
}
void polymulti(int n, int m, complex F[], complex G[])
{
    int limit, i;
    for (limit = 1; limit <= n + m; limit <<= 1, ++len)
        ;
    for (i = 0; i < limit; ++i)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (len - 1)); //二进制反转
    FFT(limit, F, 1);
    FFT(limit, G, 1);
    for (int i = 0; i < limit; ++i)
        F[i] = F[i] * G[i];
    FFT(limit, F, -1);
    for (i = 0; i <= n + m; ++i)
        F[i].x = (int)(F[i].x / limit + 0.5);
}

void divide(complex g[],complex f[],int L,int R){
    if (L==R) return;
    int mid = (L+R)>>1;
    divide(g,f,L,mid);
    for (int i = 0; i <= (R - L + 1); ++i)
        a[i] = b[i] = 0;
    for (int i=L;i<=mid;++i) a[i-L] = f[i];
    for (int i=0;i<=R-L;++i) b[i] = g[i];
    polymulti(R-L+1,R-L+1,a,b);
    for (int i=mid+1;i<=R;++i) f[i].x = md(int(f[i].x+a[i-L].x));
    divide(g,f,mid+1,R);
}

int main()
{
    scanf("%d",&n);
    int t;
    for (int i=0;i<n;++i){
        scanf("%d",&t);
        g[i].x = t;
    }
    f[0] = 1;
    divide(f,g,1,n-1);
    return 0;
}
#include <cstdio>
#include <cmath>
#include <iomanip>
#include <algorithm>

const double pi = acos(-1.0);
const int MAXN = 4 * 1e6 + 10;
int n, m, limit, i, len(0), r[MAXN * 2];
struct complex
{
    double x, y;
    complex(double xx = 0, double yy = 0) : x(xx), y(yy) {}
};
complex operator+(complex a, complex b) { return complex(a.x + b.x, a.y + b.y); }
complex operator-(complex a, complex b) { return complex(a.x - b.x, a.y - b.y); }
complex operator*(complex a, complex b) { return complex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }

complex F[MAXN], G[MAXN];

inline int readint()
{
    int ret(0), sgn(1);
    char c;
    while (isspace(c = getchar()))
        ;
    if (c == '-')
    {
        sgn = -1;
        c = getchar();
    }
    while (isdigit(c))
    {
        ret = (ret << 3) + (ret << 1) + c - '0';
        c = getchar();
    }
    return ret * sgn;
}

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

int main()
{

    n = readint();
    m = readint();
    for (i = 0; i <= n; ++i)
        F[i].x = readint();
    for (i = 0; i <= m; ++i)
        G[i].x = readint();
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
        printf("%d ", (int)(F[i].x / limit + 0.5));
    return 0;
}

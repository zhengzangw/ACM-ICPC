#ifndef math_H
#define math_H

#include <algorithm>
#include <cassert>
#include <string>
#include <cmath>
#include <map>

using int64 = long long;
#define MAXN_ 10000
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) ((a)>(b))? (b) : (a))
#define eps 1e-10
#define fabs(x) ((x) > 0 ? (x) : -(x))
#define randint(x, y) rand() % (y) + (x)
#define md(x) ((x) % mod + mod) % mod
const double pi = acos(-1.0);

int FM(int a, int b, int mod) //快速幂
{
  int c = 1;
  for (; b; b >>= 1, a = 1ll * a * a % mod)
    if (b & 1)
      c = 1ll * c * a % mod;
  return c;
}

int ni(int a, int mod)
{
  return FM(a, mod - 2, mod);
}

int FM(int a, int b)
{
  int c = 1;
  for (; b; b >>= 1, a = 1ll * a * a)
    if (b & 1)
      c = 1ll * c * a;
  return c;
}

int multiMod(int a, int b, int mod) //快速乘
{
  b %= mod;
  a %= mod;
  int res = 0;
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

int bigMod(std::string val, int mod) //长整数取余
{
  int res = 0;
  for (int i = 0; i < val.length(); ++i)
  {
    res = ((res)*10 + val[i] - '0') % mod;
  }
  return res;
}
int gcd(int a, int b)
{
  if (a < b)
    std::swap(a, b);
  while (b != 0)
  {
    int c = a;
    a = b;
    b = c % b;
  }
  return a;
}

int lcm(int a, int b)
{
  int g = gcd(a, b);
  return a / g * b;
}

int ext_gcd(int a, int b, int &x, int &y)
{
  if (b == 0)
  {
    x = 1;
    y = 0;
    return a;
  }
  else
  {
    int r = ext_gcd(b, a % b, y, x);
    y -= x * (a / b);
    return r;
  }
}

int ni_mod(int a, int mod)
{ // only for (x,mod)=1
  int x, y;
  ext_gcd(mod, a, x, y);
  return ((y % mod) + mod) % mod;
}

int prime_table(int n) //线性筛法
{
  bool p[MAXN_];
  int cprime[MAXN_];
  int used = 0;
  p[0] = 1;
  p[1] = 1;
  for (int i = 2; i <= n; ++i)
  {
    if (!p[i])
    {
      cprime[used++] = i;
      //积性函数f(p)
    }
    for (int j = 0; j < used; ++j)
    {
      if (i * cprime[j] > n)
        break;
      p[i * cprime[j]] = true;
      if (i % cprime[j] == 0)
      { // n=p1 * (p2*p3*...*pk)
        break;
        //积性函数f(p*p^k*n)=F(f(p))*f(p^k*n)
      }
      //积性函数f(mn)=f(m)f(n)
    }
  }
  return used;
}

//Based on 埃式筛法O(nlgnlgn)
void segment_prime(int a, int b, bool p[], bool sqrtbp[])
{
  for (int i = 2; i * i <= b; ++i)
    if (!sqrtbp[i])
    {
      for (int j = i * i; j * j <= b; j += i)
        sqrtbp[j] = 1;
      for (int j = max(i * i, (a - 1) / i + 1) * i; j <= b; j += i)
        p[j - a] = 1;
    }
}

int euler_phi(int n)
{
  int ans = n;
  for (int i = 2; i * i <= n; ++i)
    if (ans % i == 0)
    {
      ans = ans / i * (i - 1);
      while (n % i == 0)
        n /= i;
    }
  if (n > 1)
    ans = ans / n * (n - 1);
  return ans;
}

void euler_table(int n, int euler[]) //O(nlgnlgn)
{
  euler[1] = 1;
  for (int i = 2; i < n; i++)
    euler[i] = i;
  for (int i = 2; i < n; i++)
    if (euler[i] == i)
      for (int j = i; j < n; j += i)
        euler[j] = euler[j] / i * (i - 1);
}

bool cong_eq(int a, int b, int m, int ans[], int &len)
{
  int d, x, y;
  d = ext_gcd(a, m, x, y);
  if (b % d)
    return false;
  int base = ((b / d * x) % m + m) % m;
  len = d;
  for (int i = 0; i < len; ++i)
    ans[i] = (base + i * (m / len)) % m;
  return true;
}

int BSGSx(int A, int C, int mod) //离散对数
{
  A %= mod;
  C %= mod;
  if (C == 1)
    return 0;
  int cnt = 0, tmp = 1;

  for (int g = gcd(A, mod); g != 1; g = gcd(A, mod))
  {
    if (C % g)
      return -1;
    C /= g;
    mod /= g;
    tmp = tmp * A / g % mod;
    ++cnt;
    if (C == tmp)
      return cnt;
  }
  //tmp*a^(x-cnt)=b' (mod c')

  int T = (int)sqrt(0.5 + mod);
  int b = C;
  std::map<int, int> hash;
  hash[b] = 0;
  for (int i = 1; i <= T; ++i)
  {
    b = b * A % mod;
    hash[b] = i;
  }

  A = FM(A, T, mod);
  for (int u = 1; u <= T; ++u)
  {
    tmp = tmp * A % mod;
    if (hash.count(tmp))
      return u * T - hash[tmp] + cnt;
  }
  return -1;
}

int china(int a[], int b[], int r)
{
  // x = a[i] mod b[i];
  int M = 1, i, Mi, ans = 0;
  for (i = 0; i < r; ++i)
    M *= b[i];
  for (i = 0; i < r; ++i)
  {
    Mi = M / b[i];
    ans = (ans + ni_mod(Mi, b[i]) * Mi * a[i]) % M;
  }
  return (ans + M) % M;
}

int cong_eq_group(int a[], int b[], int r)
{
  int ans = 0, M = 1, x, y, d;
  for (int i = 0; i < r; ++i)
  {
    d = ext_gcd(M, b[i], x, y); //b[i]*y=d(mod M)
    if ((ans - a[i]) % d)
      return -1;
    M = M / d * b[i];
    int z = y * ((ans - a[i]) / d);
    ans = z * b[i] + a[i];
    ans = ((ans % M) + M) % M;
  }
  return ans;
}

bool witness(int s, int n)
{
  int u = n - 1;
  int t = 0;
  while ((u & 1) == 0)
    u >>= 1, t++; //二次探测定理+Fermat测试

  int x = FM(s, u, n), tmp;
  while (t--)
  {
    tmp = x;
    x = multiMod(x, x, n);
    if (x == 1)
    {
      if (tmp == n - 1 || tmp == 1)
        return false;
      else
        return true;
    }
  }
  return true;
}

bool millerRabin(int n, const int times = 3)
{ //P=4^(-times)
  if (n == 2)
    return true;
  if ((n & 1) == 0 || n < 2)
    return false;
  int i = times;
  while (i--)
  {
    int s = randint(1, n - 1);
    if (witness(s, n))
      return false; //s^(n-1)==1 (mod n)
  }
  return true;
}

int pollard_rho(int n)
{ //O(n^(1/4)) -- Why?
  int x, y, k = 2, d, i = 1, c;
  x = y = randint(0, n - 1);
  c = randint(0, n - 1);
  while (true)
  {
    ++i;
    x = (multiMod(x, x, n) + c) % n;
    if (y == x)
      return 1;
    else if (y > x)
      d = gcd(y - x, n);
    else
      d = gcd(x - y, n);
    if (d != 1 && d != n - 1)
      return d;
    else
    {
      if (i == k)
      { // Floyd’s cycle detection trick
        y = x;
        k <<= 1;
      }
    }
  }
}

void split(int n, int factors[], int &len)
{
  if (millerRabin(n))
    factors[++len] = n; //or ++fac[n] or (map<int,int> ++fac[n])
  else
  {
    int p;
    do
    {
      p = pollard_rho(n);
    } while (p == n || p == 1);
    split(p, factors, len);
    split(n / p, factors,len);
  }
}

//O(sqrt(n)/logn)
int prime[1];
int split_(int n, int a[], int &len)
{
  int sqrtn = (int)sqrt(n);
  len = 0;
  a[0] = 1;
  int ans = 1, tmp;
  for (int i = 0; prime[i] <= sqrtn; ++i)
  {
    tmp = 0;
    while (n % prime[i] == 0)
    {
      tmp++;
      a[++len] = prime[i];
      n /= prime[i];
    }
    ans *= tmp + 1;
  }
  if (n != 1)
  {
    a[++len] = n;
    ans <<= 1;
  }
  return ans;
}

void split__(int n, int a[], int &len)
{
  int sqrtn = (int)sqrt(n);
  len = 0;
  a[0] = 1;
  if ((n & 1) == 0)
  {
    a[++len] = 2;
    n /= 2;
  }
  for (int i = 3; i <= sqrtn; i += 2)
    if (n % i == 0)
    {
      a[++len] = i;
      while (n % i == 0)
        n /= i;
    }
  if (n > 1)
    a[++len] = n;
}

int ord(int a, int m)
{
  int phi = euler_phi(m), tmp = 1;
  for (int i = 1; i <= phi; i++)
  {
    tmp = tmp * a % m;
    if (tmp == 1)
      return i;
  }
  return -1;
}

int is_g(int x, int factors[], int &len, int phi)
{
  for (int i = 0; i < len; ++i)
    if (factors[i] != 1 && FM(x, phi / factors[i], phi) == 1)
      return false;
  return true;
}

int get_g(int p)
{
  if (p == 2 || p == 4)
    return 1;
  if (!(millerRabin(p) || (!(p & 1) && millerRabin(p / 2))))
    return -1;
  int phi = euler_phi(p), factors[1000], len = 0;
  split(phi, factors, len);
  for (int i = 1; i < p; i++)
    if (is_g(i, factors, len, phi))
      return i;
  return -1;
}

template <class T>
struct Matrix
{
  T m[MAXN_][MAXN_];
  int l, r;
  Matrix<T>(int w)
  {
    l = r = w;
    memset(m, 0, sizeof(m));
  }
  Matrix<T>()
  {
    l = r = 0;
    memset(m, 0, sizeof(m));
  }

  friend Matrix operator*(const Matrix &a, const Matrix &b)
  {
    Matrix c;
    memset(c.m, 0, sizeof c.m);
    c.l = a.l, c.r = b.r;
    for (int i = 0; i < a.l; ++i)
      for (int j = 0; j < b.r; ++j)
        for (int k = 0; k < a.r; ++k)
          c.m[i][j] = (c.m[i][j] + a.m[i][k] * a.m[k][j]);
    return c;
  }

  friend Matrix operator*=(Matrix &a, const Matrix &b)
  {
    a = a * b;
    return a;
  }

  friend Matrix operator+(const Matrix &a, const Matrix &b)
  {
    assert(a.l == b.l);
    assert(a.r == b.r);
    Matrix c;
    c.l = a.l;
    c.r = a.r;
    for (int i = 0; i < c.l; ++i)
      for (int j = 0; j < c.r; ++j)
        c.m[i][j] = a.m[i][j] + b.m[i][j];
    return c;
  }

  friend Matrix operator+=(Matrix &a, const Matrix &b)
  {
    a = a + b;
    return a;
  }

  friend Matrix operator^(const Matrix a, int t)
  {
    Matrix ans;
    ans.l = ans.r = a.l;
    for (int i = 0; i < a.l; i++)
      ans.m[i][i] = 1;
    Matrix tmp = a;
    while (t)
    {
      if (t & 1)
        ans = tmp * a;
      tmp = tmp * tmp;
      t >>= 1;
    }
    return ans;
  }
};

char gauss_cpivot(int n, double a[][MAXN_], double b[])
{
  // Solve A*X=b, Ans is in b[] O(n^3)
  int i, j, k, row;
  double maxp, t;
  // For every column
  for (k = 0; k < n; k++)
  {
    // Find the max element
    for (maxp = 0, i = k; i < n; i++)
      if (fabs(a[i][k]) > fabs(maxp))
        maxp = a[row = i][k];
    // If failed, then there is no only solution
    if (fabs(maxp) < eps)
      return 0;
    // Exchange Row row and Row k
    if (row != k)
    {
      for (j = k; j < n; j++)
        t = a[k][j], a[k][j] = a[row][j], a[row][j] = t;
      t = b[k], b[k] = b[row], b[row] = t;
    }
    // Update rest elements
    for (j = k + 1; j < n; j++)
    {
      a[k][j] /= maxp;
      for (i = k + 1; i < n; i++)
        a[i][j] -= a[i][k] * a[k][j];
    }
    b[k] /= maxp;
    for (i = k + 1; i < n; i++)
      b[i] -= b[k] * a[i][k];
  }
  // Deal with the 上三角矩阵
  for (i = n - 1; i >= 0; i--)
    for (j = i + 1; j < n; j++)
      b[i] -= a[i][j] * b[j];
  return 1;
}

int valid_multi(int a, int b)
{
  int MAX = (int)1e18;
  //两数相乘是否越界
  int ret = 0;
  for (; b; b >>= 1)
  {
    if (b & 1)
    {
      ret += a;
      if (ret > MAX || a < 0)
        return 0;
    }
    a += a;
    if (a >= MAX)
      a = -1;
  }
  return 1;
}

//欧拉降幂公式
int FM_Euler(int a, int b, int m)
{
  int p = euler_phi(m);
  b = b % p;
  return FM(a, b, m);
}

//复数
struct complex
{
  double x, y;
  complex(double xx = 0, double yy = 0) : x(xx), y(yy) {}
};
complex operator+(complex a, complex b) { return complex(a.x + b.x, a.y + b.y); }
complex operator-(complex a, complex b) { return complex(a.x - b.x, a.y - b.y); }
complex operator*(complex a, complex b) { return complex(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x); }
//FFT O(nlogn) 多项式乘法，两个数组卷积
//分治法慢 A(w)=A'(w^2)+wA(w^2) A(w+n/2) = A'(w^2)-wA(w^2)
//迭代实现 10e6
//卷积 ck = \sum_{i+j=k}(ai+aj)
int r[MAXN_ * 2];
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
typedef long long LL;
const int P = 998244353, G = 3, Gni = 332748118;
//998244353 = 119*2^23 + 1; g = 3;
//1004535809 = 479*2^21 + 1; g = 3;
//2281701377 = 17*2^27 + 1; g = 3;
void NTT(int limit, LL *A, int type)
{
  for (int i = 0; i < limit; ++i)
    if (i < r[i])
      std::swap(A[i], A[r[i]]);
  for (int mid = 1; mid < limit; mid <<= 1)
  {
    LL wn = FM(type == 1 ? G : Gni, (P - 1) / (mid << 1), P);
    for (int j = 0; j < limit; j += (mid << 1))
    {
      LL w = 1;
      for (int k = 0; k < mid; k++, w = (w * wn) % P)
      {
        LL x = A[j + k], y = w * A[j + k + mid] % P;
        A[j + k] = (x + y) % P;
        A[j + k + mid] = (x - y + P) % P;
      }
    }
  }
}
void polymulti(int n, int m, complex F[], complex G[])
{
  int limit, i,len;
  for (limit = 1,len = 0; limit <= n + m; limit <<= 1, ++len)
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
void polymulti(int n,int m,LL F[],LL G[])
{
  int limit,i,len;
  for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
    ;
  for (int i = 0; i < limit; ++i)
    r[i] = (r[i >> 1] >> 1 | ((i & 1) << (len - 1)));
  NTT(limit, F, 1);
  NTT(limit, G, 1);
  for (int i = 0; i < limit; i++)
    F[i] = (F[i] * G[i]) % P;
  NTT(limit,F, -1);
  LL inv = FM(limit, P - 2, P);
  for (int i = 0; i <= n + m; ++i)
    F[i] = (F[i] * inv) % P;
}

//求1-n中与a互质的数个数
//容斥原理
int prime_of_n[100];
int mutual_prime_(int a, int n)
{
  int len;
  split(a, prime_of_n, len);
  int mask = 1 << len, ans = 0;
  for (int i = 1; i < mask; ++i)
  {
    int cnt = 0, v = 1;
    for (int j = 0; j < len; ++j)
      if (i & (1 << j))
      {
        cnt++;
        v *= prime_of_n[j + 1];
      }
    if (cnt & 1)
      ans += n / v;
    else
      ans -= n / v;
  }
  return n - ans;
}
//莫比乌斯反演
int mu[100];
int mutual_prime(int a,int b)
{
  int ans;
  if (a>b) std::swap(a,b);
  for (int i = 1; i <= a; ++i)
    ans += mu[i] * (a / i) * (b / i);
}

#endif
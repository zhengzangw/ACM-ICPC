#ifndef math_H
#define math_H

#include <algorithm>
#include <cassert>
#include <cmath>
#include <map>
#include <string>

typedef long long LL;
#define MAXN 10000
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) (((a)>(b))? (b) : (a))
#define eps 1e-10
#define fabs(x) ((x) > 0 ? (x) : -(x))
#define randint(x, y) rand() % (y) + (x)
#define md(x) (((1ll * x) % mod + mod) % mod)
const int mod = 998244353;
const double pi = acos(-1.0);

//快速幂
int FM(int a, int b, int mod) {
    int c = 1;
    for (; b; b >>= 1, a = 1ll * a * a % mod)
        if (b & 1) c = 1ll * c * a % mod;
    return c % mod;
}

//快乘
int multiMod(int a, int b, int mod) {
    b %= mod;
    a %= mod;
    int res = 0;
    while (b) {
        if (b & 1) {
            res += a;
            if (res >= mod) res -= mod;
        }
        b >>= 1;
        a = a << 1;
        if (a >= mod) a -= mod;
    }
    return res;
}

//长整数取余
int bigMod(std::string val, int mod) {
    int res = 0;
    for (int i = 0; i < val.length(); ++i) {
        res = ((res)*10 + val[i] - '0') % mod;
    }
    return res;
}

//欧几里得算法
int gcd(int a, int b) {
    if (a < b) std::swap(a, b);
    while (b != 0) {
        int c = a;
        a = b;
        b = c % b;
    }
    return a;
}

//最小公倍数
int lcm(int a, int b) {
    int g = gcd(a, b);
    return a / g * b;
}

//扩展欧几里得
int ext_gcd(int a, int b, int& x, int& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    } else {
        int r = ext_gcd(b, a % b, y, x);
        y -= x * (a / b);
        return r;
    }
}
//逆元
//扩展欧几里得
int inv(int a, int mod) {  // only for (x,mod)=1
    a %= mod;
    int x, y;
    ext_gcd(mod, a, x, y);
    return md(y);
}
// FM(a,mod-2,mod) mod为素数 O(log(mod))
//预处理fac,fac_inv逆元
int fac[MAXN], fac_inv[MAXN],n_inv[MAXN];
void pre_fac_inv(int n) {
    fac[0] = 1;
    for (int i = 1; i <= n; ++i) {
        fac[i] = md(1ll * fac[i - 1] * i);
    }
    fac_inv[n] = inv(fac[n], mod);
    for (int i = n - 1; i >= 0; ++i) {
        fac[i] = md(1ll * fac[i + 1] * (i + 1));
    }
}
//预处理前n个数逆元
void pre_n_inv(int n) {
    n_inv[1] = 1;
    for (int i = 2; i <= n; i++) {
        n_inv[i] = md(1ll * md(-mod / i) * md(n_inv[mod % i]));
    }
    /* 更新fac_inv另一种方法
    fac_inv[0] = fac_inv[1] = 1;
    for (int i= 2;i<=n;i++){
        fac_inv[i] = fac_inv[i-1]*n_inv[i];
    }*/
}

//lucas 算法
//O(log_p n) C_n^m
LL lucas(LL n,LL m,LL mod){
    //mod is prime
    if (n<m) return 0;
    else if (n<mod) return md(md(fac[n]*fac_inv[m])*fac_inv[n-m]);
    else return md(lucas(n/mod,m/mod,mod)*lucas(n%mod,m%mod,mod));
}

//线性筛法
int prime_table(int n) {
    bool vis[MAXN];
    int prime[MAXN];
    int used = 0;
    vis[0] = vis[1] = 1;
    for (int i = 2; i <= n; ++i) {
        if (!vis[i]) {
            prime[used++] = i;
            //积性函数f(p)
            // phi(p) = p-1
        }
        for (int j = 0; j < used; ++j) {
            if (i * prime[j] > n) break;
            vis[i * prime[j]] = true;
            if (i % prime[j] == 0) {  // n=p1 * (p2*p3*...*pk)
                //积性函数f(p*p^k*n)=F(f(p))*f(p^k*n)
                // phi(i*prime[j]) = prime[j]*phi(i)
                break;
            }
            //积性函数f(mn)=f(m)f(n)
            // phi(i*prime[j]) = (prime[j]-1)*phi(i)
        }
    }
    return used;
}

// Based on 埃式筛法O(nlgnlgn)
void segment_prime(int a, int b, bool p[], bool sqrtbp[]) {
    for (int i = 2; i * i <= b; ++i)
        if (!sqrtbp[i]) {
            for (int j = i * i; j * j <= b; j += i) sqrtbp[j] = 1;
            for (int j = max(i * i, (a - 1) / i + 1) * i; j <= b; j += i)
                p[j - a] = 1;
        }
}

//单点欧拉函数
int phi(int n) {
    int ans = n;
    for (int i = 2; i * i <= n; ++i)
        if (ans % i == 0) {
            ans = ans / i * (i - 1);
            while (n % i == 0) n /= i;
        }
    if (n > 1) ans = ans / n * (n - 1);
    return ans;
}

//欧拉降幂公式
int FM_Euler(int a, int b, int m) {
    int p = phi(m);
    b = b % p;
    return FM(a, b, m);
}

//同余方程
bool cong_eq(int a, int b, int m, int ans[], int& len) {
    int d, x, y;
    d = ext_gcd(a, m, x, y);
    if (b % d) return false;
    int base = ((b / d * x) % m + m) % m;
    len = d;
    for (int i = 0; i < len; ++i) ans[i] = (base + i * (m / len)) % m;
    return true;
}

//同余方程组 中国剩余定理
int china(int a[], int b[], int r) {
    // x = a[i] mod b[i];
    int M = 1, i, Mi, ans = 0;
    for (i = 0; i < r; ++i) M *= b[i];
    for (i = 0; i < r; ++i) {
        Mi = M / b[i];
        ans = (ans + inv(Mi, b[i]) * Mi * a[i]) % M;
    }
    return (ans + M) % M;
}

//一般同余方程组
int cong_eq_group(int a[], int b[], int r) {
    int ans = 0, M = 1, x, y, d;
    for (int i = 0; i < r; ++i) {
        d = ext_gcd(M, b[i], x, y);  // b[i]*y=d(mod M)
        if ((ans - a[i]) % d) return -1;
        M = M / d * b[i];
        int z = y * ((ans - a[i]) / d);
        ans = z * b[i] + a[i];
        ans = ((ans % M) + M) % M;
    }
    return ans;
}

//离散对数
int BSGSx(int A, int C, int mod) {
    A %= mod;
    C %= mod;
    if (C == 1) return 0;
    int cnt = 0, tmp = 1;

    for (int g = gcd(A, mod); g != 1; g = gcd(A, mod)) {
        if (C % g) return -1;
        C /= g;
        mod /= g;
        tmp = tmp * A / g % mod;
        ++cnt;
        if (C == tmp) return cnt;
    }
    // tmp*a^(x-cnt)=b' (mod c')

    int T = (int)sqrt(0.5 + mod);
    int b = C;
    std::map<int, int> hash;
    hash[b] = 0;
    for (int i = 1; i <= T; ++i) {
        b = b * A % mod;
        hash[b] = i;
    }

    A = FM(A, T, mod);
    for (int u = 1; u <= T; ++u) {
        tmp = tmp * A % mod;
        if (hash.count(tmp)) return u * T - hash[tmp] + cnt;
    }
    return -1;
}

// MillerRabin 素数测试
bool witness(int s, int n) {
    int u = n - 1;
    int t = 0;
    while ((u & 1) == 0) u >>= 1, t++;  //二次探测定理+Fermat测试

    int x = FM(s, u, n), tmp;
    while (t--) {
        tmp = x;
        x = multiMod(x, x, n);
        if (x == 1) {
            if (tmp == n - 1 || tmp == 1)
                return false;
            else
                return true;
        }
    }
    return true;
}

bool millerRabin(int n, const int times = 3) {  // P=4^(-times)
    if (n == 2) return true;
    if ((n & 1) == 0 || n < 2) return false;
    int i = times;
    while (i--) {
        int s = randint(1, n - 1);
        if (witness(s, n)) return false;  // s^(n-1)==1 (mod n)
    }
    return true;
}

// pollard 质因数分解
int pollard_rho(int n) {  // O(n^(1/4)) -- Why?
    int x, y, k = 2, d, i = 1, c;
    x = y = randint(0, n - 1);
    c = randint(0, n - 1);
    while (true) {
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
        else {
            if (i == k) {  // Floyd’s cycle detection trick
                y = x;
                k <<= 1;
            }
        }
    }
}

void split(int n, int factors[], int& len) {
    if (millerRabin(n))
        factors[++len] = n;  // or ++fac[n] or (map<int,int> ++fac[n])
    else {
        int p;
        do {
            p = pollard_rho(n);
        } while (p == n || p == 1);
        split(p, factors, len);
        split(n / p, factors, len);
    }
}

//普通质因数分解 O(sqrt(n)/logn)
//统计因子个数
int prime[1];
int split_(int n, int a[], int& len) {
    int sqrtn = (int)sqrt(n);
    len = 0;
    a[0] = 1;
    int ans = 1, tmp;
    for (int i = 0; prime[i] <= sqrtn; ++i) {
        tmp = 0;
        while (n % prime[i] == 0) {
            tmp++;
            a[++len] = prime[i];
            n /= prime[i];
        }
        ans *= tmp + 1;
    }
    if (n != 1) {
        a[++len] = n;
        ans <<= 1;
    }
    return ans;
}

//阶
int ord(int a, int m) {
    int p = phi(m), tmp = 1;
    for (int i = 1; i <= p; i++) {
        tmp = tmp * a % m;
        if (tmp == 1) return i;
    }
    return -1;
}

//原根
int is_g(int x, int factors[], int& len, int phi) {
    for (int i = 0; i < len; ++i)
        if (factors[i] != 1 && FM(x, phi / factors[i], phi) == 1) return false;
    return true;
}

int get_g(int p) {
    if (p == 2 || p == 4) return 1;
    if (!(millerRabin(p) || (!(p & 1) && millerRabin(p / 2)))) return -1;
    int phiofn = phi(p), factors[1000], len = 0;
    split(phiofn, factors, len);
    for (int i = 1; i < p; i++)
        if (is_g(i, factors, len, phiofn)) return i;
    return -1;
}

//矩阵
template <class T>
struct Matrix {
    T m[MAXN][MAXN];
    int l, r;
    Matrix<T>(int w) {
        l = r = w;
        memset(m, 0, sizeof(m));
    }
    Matrix<T>() {
        l = r = 0;
        memset(m, 0, sizeof(m));
    }

    void print(){
        for (int i=0;i<l;++i,puts("")){
            printf("%d",m[i][0]);
            for (int j=1;j<r;++j)
                printf(" %d",m[i][j]);
        }
    }

    friend Matrix operator*(const Matrix& a, const Matrix& b) {
        Matrix c;
        memset(c.m, 0, sizeof c.m);
        c.l = a.l, c.r = b.r;
        for (int i = 0; i < a.l; ++i)
            for (int j = 0; j < b.r; ++j)
                for (int k = 0; k < a.r; ++k)
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
        for (int i = 0; i < c.l; ++i)
            for (int j = 0; j < c.r; ++j) c.m[i][j] = a.m[i][j] + b.m[i][j];
        return c;
    }

    friend Matrix operator+=(Matrix& a, const Matrix& b) {
        a = a + b;
        return a;
    }

    friend Matrix operator^(const Matrix a, int t) {
        Matrix ans;
        ans.l = ans.r = a.l;
        for (int i = 0; i < a.l; i++) ans.m[i][i] = 1;
        Matrix tmp = a;
        while (t) {
            if (t & 1) ans = tmp * a;
            tmp = tmp * tmp;
            t >>= 1;
        }
        return ans;
    }
};

//矩阵求逆 O(n^3)
template <typename T>
Matrix<T> inv(Matrix<T> a,int n){
    int m = n<<1;
    for (int i=0;i<n;i++)
        a[i][n+i] = 1;
    for (int k=0;k<n;k++){
        for (int i = k;i<n;i++)
            if (a[i][k]){
                for (int j = 1;j<m;j++) swap(a[k][j],a[i][j]);
                break;
            }
        if (!a[k][k]) {puts("No Solution"); return 0;}
        LL r = FM(a[k][k],mod-2,mod);
        for (int j=k;j<m;++j) a[k][j] = md(a[k][j]*r);
        for (int i=0;i<n;i++){
            r = a[i][k];
            if (i!=k)
            for (int j = k; j < m; j++) a[i][j] = md(a[i][j] - md(r*a[k][j]));
        }
    }
    Matrix<T> a_inv(n);
    for (int i=0;i<n;i++)
        for (int j=0;j<n;j++)
            a_inv[i][j] = a[i][n+j];
    return a_inv;
}

//高斯消元
bool gauss_cpivot(int n, double a[][MAXN], double b[]) {
    // Solve A*X=b, Ans is in b[] O(n^3)
    int i, j, k, row;
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
            for (j = k; j < n; j++){
                t = a[k][j]; a[k][j] = a[row][j]; a[row][j] = t;
            }
            t = b[k]; b[k] = b[row]; b[row] = t;
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

//判断两数相乘是否越界
int valid_multi(int a, int b) {
    int MAX = (int)1e18, ret = 0;
    for (; b; b >>= 1) {
        if (b & 1) {
            ret += a;
            if (ret > MAX || a < 0) return 0;
        }
        a += a;
        if (a >= MAX) a = -1;
    }
    return 1;
}

//复数
struct complex_ {
    double x, y;
    complex_(double xx = 0, double yy = 0) : x(xx), y(yy) {}
    friend complex_ operator+(complex_ a, complex_ b) {
        return complex_(a.x + b.x, a.y + b.y);
    }
    friend complex_ operator-(complex_ a, complex_ b) {
        return complex_(a.x - b.x, a.y - b.y);
    }
    friend complex_ operator*(complex_ a, complex_ b) {
        return complex_(a.x * b.x - a.y * b.y, a.x * b.y + a.y * b.x);
    }
};
// typedef complex<double> complex_

// 998244353 = 119*2^23 + 1; g = 3;gni2 = 332748118;
// 1004535809 = 479*2^21 + 1; g = 3; gni3 = 334845270;
// 2281701377 = 17*2^27 + 1; g = 3;
// 469762049 = 7*2^26 + 1; g = 3; gni1 = 156587350;


//NTT O(nlogn)
//For MTT
/*
int Pnum = 1;
LL Proot_cur,Proot_invcur,Pcur;
void set_prime(){
    Proot_cur=Proot[Pnum];
    Proot_invcur=Proot_inv[Pnum];
    Pcur=P[Pnum];
}
*/
int r[MAXN * 2];
void calrev(int n, int len) {
    for (int i = 0; i < n; ++i)
        r[i] = (r[i >> 1] >> 1) | ((i & 1) << (len - 1));
}
const LL P[] = {469762049,998244353,1004535809};
const LL Proot[] = {3,3,3};
const LL Proot_inv[] = {156587350,332748118,334845270};
const int Pnum = 1;
const LL Proot_cur=Proot[Pnum],Proot_invcur=Proot_inv[Pnum],Pcur=P[Pnum];

void NTT(int limit,LL* A, int type) {
    for (int i = 0; i < limit; ++i)
        if (i < r[i]) std::swap(A[i], A[r[i]]);
    for (int mid = 1; mid < limit; mid <<= 1) {
        LL wn = FM(type == 1 ? Proot_cur : Proot_invcur, (Pcur - 1) / (mid << 1), Pcur);
        for (int j = 0; j < limit; j += (mid << 1)) {
            LL w = 1;
            for (int k = 0; k < mid; k++, w = (w * wn) % Pcur) {
                LL x = A[j + k], y = w * A[j + k + mid] % Pcur;
                A[j + k] = (x + y) % Pcur;
                A[j + k + mid] = (x - y + Pcur) % Pcur;
            }
        }
    }
    if (type==-1){
        LL inv = FM(limit, Pcur - 2, Pcur);
        for (int i = 0; i < limit; ++i) A[i] = (A[i] * inv) % Pcur;
    }
}

LL sa[MAXN],sb[MAXN];
void poly_mult(int n, int m, LL a[], LL b[],LL ans[]) {
    memset(sa,0,sizeof(sa)); memset(sb,0,sizeof(sb));
    //set_prime();
    int limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    for (int i=0;i<limit;++i) {sa[i] = a[i]; sb[i] = b[i];}
    calrev(limit, len);
    NTT(limit, sa, 1);
    NTT(limit, sb, 1);
    for (int i = 0; i < limit; i++) ans[i] = (sa[i] * sb[i]) % Pcur;
    NTT(limit, ans, -1);
}

//多项式求逆 mod x^n
LL t1[MAXN], t2[MAXN];
void poly_inv(LL* a, LL n, LL mod,LL* ans) {
    memset(t1,0,sizeof(t1)); memset(t2,0,sizeof(t2));
    //set_prime();
    ans[0] = FM(a[0], mod - 2, mod);
    for (int bas = 2, len = 2, limit = 4; bas < (n << 1); bas <<= 1, ++len, limit<<=1) {
        calrev(limit, len);
        for (int i = 0; i < bas; ++i) {
            t1[i] = a[i]; t2[i] = ans[i];
        }
        NTT(limit, t1, 1);
        NTT(limit, t2, 1);
        for (int i = 0; i < limit; ++i)
            ans[i] = md(t2[i] * md(2 - md(t1[i] * t2[i])));
        NTT(limit,ans,-1);
        for (int i = 0; i < bas; ++i) ans[i] %= mod;
        for (int i=bas; i <limit;++i) ans[i] = 0;
    }
}

//多项式除法
//O(nlogn) F(x) = Q(x)G(x) + R(x) n = m*(n-m) + (m-1)
LL G_R[MAXN],G_R_inv[MAXN],F_R[MAXN],Q_R[MAXN],G_Q[MAXN];
void poly_div(int n,int m,int mod,LL F[],LL G[],LL Q[],LL R[]){
    for (int i=0;i<=m;++i) G_R[i] = G[m-i];
    for (int i=n-m+1;i<=n;++i) G_R[i] = 0;
    poly_inv(G_R, n-m+1, mod, G_R_inv);
    for (int i=0;i<=n-m;++i) F_R[i] = F[n-i];
    poly_mult(n-m, n-m, F_R, G_R_inv, Q_R);
    for (int i=0;i<=n-m;++i) Q[i] = Q_R[n-m-i];
    poly_mult(m,n-m,G,Q,G_Q);
    for (int i=0;i<=n;++i) R[i] = md(F[i] - G_Q[i]);
}

//任意模数NTT
const LL P1 = P[0],P2 = P[1],P3 = P[2];
const LL M = P1 * P2;
LL a1[MAXN], a2[MAXN], a3[MAXN], b1[MAXN], b2[MAXN], b3[MAXN];
void MTT(int n, int m, LL mod, LL a[], LL b[], LL ans[]) {
    for (int i = 0; i <= n; ++i) a1[i] = a2[i] = a3[i] = a[i];
    for (int i = 0; i <= m; ++i) b1[i] = b2[i] = b3[i] = b[i];
    int limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    calrev(limit, len);
    //Pnum = 0; 
    poly_mult(n, m, a1, b1, a1);
    //Pnum = 1; 
    poly_mult(n, m, a2, b2, a2);
    //Pnum = 2; 
    poly_mult(n, m, a3, b3, a3);
    for (int i = 0; i <= n + m; ++i) {
        LL tmp = (multiMod((P2 * a1[i]), inv(P2, P1), M) +
                  multiMod((P1 * a2[i]) % M, inv(P1, P2), M)) % M;
        LL k = ((a3[i] - tmp) % P3 + P3) % P3 * inv(M, P3) % P3;
        ans[i] = ((k % mod) * (M % mod) % mod + tmp % mod) % mod;
    }
}

// FFT O(nlogn) 多项式乘法，两个数组卷积
//分治法慢 A(w)=A'(w^2)+wA(w^2) A(w+n/2) = A'(w^2)-wA(w^2)
//迭代实现 \e6
//卷积 ck = \sum_{i+j=k}(ai+aj)
void FFT(int limit, complex_* a, int type) {
    complex_ x, y;
    for (int i = 0; i < limit; ++i)
        if (i < r[i]) std::swap(a[i], a[r[i]]);
    for (int mid = 1; mid < limit; mid <<= 1) {
        complex_ Wn(cos(pi / mid), type * sin(pi / mid)), x, y;
        for (int j = 0; j < limit; j += (mid << 1)) {
            complex_ w(1, 0);
            for (int k = 0; k < mid; ++k, w = w * Wn) {
                x = a[j + k];
                y = w * a[j + mid + k];
                a[j + k] = x + y;
                a[j + mid + k] = x - y;
            }
        }
    }
}

complex_ ta[MAXN],tb[MAXN],tc[MAXN];
void poly_mult_double(int n, int m, int F[], int G[],int ans[]) {
    int limit, len;
    for (limit = 1, len = 0; limit <= n + m; limit <<= 1, ++len)
        ;
    for (int i=0;i<=limit;++i) {
        ta[i] = F[i]; tb[i] = G[i]; tc[i] = 0;
    }
    //二进制反转
    calrev(limit, len);
    FFT(limit, ta, 1);
    FFT(limit, tb, 1);
    for (int i = 0; i < limit; ++i) tc[i] = ta[i] * tb[i];
    FFT(limit, tc, -1);
    for (int i = 0; i <= n + m; ++i) ans[i] = (int)(tc[i].x / limit + 0.5);
}

LL TG[MAXN],TR[MAXN];
void poly_mult_mod(int m, LL *A, LL *B,LL *MOD,int opt)
{
    opt = min(opt,m<<1);
    memset(TG,0,sizeof(TG)); memset(TR,0,sizeof(TR));
    poly_mult(m,m,A,B,A);
    if (opt>m){
        poly_div(opt,m,mod,A,MOD,TG,TR);
        for (int i=m;i<=m<<1;++i) A[i] = 0;
        for (int i=0;i<m;++i) A[i] = TR[i];
    }
}


//A^n % MOD
void poly_FM(LL *A,int n,LL *ans,LL *MOD,int m){
    for (int i=1;n;n>>=1,poly_mult_mod(m,A,A,MOD,2<<i),++i){
        if (n&1) poly_mult_mod(m,ans,A,MOD,m<<1);
        if((1<<(i-1))>=(m<<1))--i;
    }
}

LL T[MAXN],C[MAXN],G[MAXN];
LL Solve_Linear_Recur(int n,int m, LL *f, LL *A){
    for (int i=1;i<=m;++i) G[m-i] = (mod-f[i])%mod;
    G[m] = 1;
    T[1] = 1; C[0] = 1;
    poly_FM(T,n,C,G,m);
    LL ans = 0;
    for (int i=0;i<=m-1;++i) ans = (ans + md(A[i] * C[i]))%mod;
    return ans;
}

//快速沃尔什变换
inline int div2(int x) { return x & 1 ? (mod + x) / 2 : x / 2; }
void FWTxor(int limit, int* a, int type) {
    int t;
    for (int mid = 1; mid < limit; mid <<= 1)
        for (int j = 0; j < limit; j += (mid << 1))
            for (int k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] = type == 1 ? md(t + a[j + k + mid])
                                     : div2(md(t + a[j + k + mid]));
                a[j + k + mid] = type == 1 ? md(t - a[j + k + mid])
                                           : div2(md(t - a[j + k + mid]));
            }
}
void FWTand(int limit, int* a, int type) {
    int t;
    for (int mid = 1; mid < limit; mid <<= 1)
        for (int j = 0; j < limit; j += (mid << 1))
            for (int k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] =
                    type == 1 ? md(t + a[j + k + mid]) : md(t - a[j + k + mid]);
                a[j + k + mid] = md(a[j + k + mid]);
            }
}
void FWTor(int limit, int* a, int type) {
    int t;
    for (int mid = 1; mid < limit; mid <<= 1)
        for (int j = 0; j < limit; j += (mid << 1))
            for (int k = 0; k < mid; ++k) {
                t = a[j + k];
                a[j + k] = md(t);
                a[j + k + mid] =
                    type == 1 ? md(t + a[j + k + mid]) : md(a[j + k + mid] - t);
            }
}

//快速莫比乌斯变换
void FMTor(int n, int* a, int type) {
    // n = log(limit)
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < 1 << n; j++)
            if (j >> i & 1)
                a[j] = md(a[j] + a[j ^ (1 << i)] * (type == 1 ? 1 : -1));
}
void FMTand(int n, int* a, int type) {
    for (int i = 0; i < n; ++i)
        for (int j = (1 << n) - 1; ~j; j--)
            if (~j >> i & 1)
                a[j] = md(a[j] + a[j | 1 << i] * (type == 1 ? 1 : -1));
}

//求1-n中与a互质的数个数
//法一：容斥原理
int prime_of_n[100];
int mutual_prime_(int a, int n) {
    int len;
    split(a, prime_of_n, len);
    int mask = 1 << len, ans = 0;
    for (int i = 1; i < mask; ++i) {
        int cnt = 0, v = 1;
        for (int j = 0; j < len; ++j)
            if (i & (1 << j)) {
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
//法二:莫比乌斯反演
int mu[100];
int mutual_prime(int a, int b) {
    int ans;
    if (a > b) std::swap(a, b);
    for (int i = 1; i <= a; ++i) ans += mu[i] * (a / i) * (b / i);
    return ans;
}

#endif

#ifndef math_H
#define math_H

#include <bits/stdc++.h>
using namespace std;
typedef long long LL;
typedef LL MATHTYPE;
#define MAXN 0
#define eps 1e-10
#define randint(x, y) rand() % (y) + (x)
#define md(x) (((1ll * x) % mod + mod) % mod)
const MATHTYPE mod = 998244353;
const double pi = acos(-1.0);


/**
 * Fast Power with mod
 * @return a^b%mod
 * @timecomplex O(nlgn)
 */
#define fp fast_power
MATHTYPE fast_power(MATHTYPE a, MATHTYPE b, MATHTYPE mod) {
    MATHTYPE c = 1;
    for (; b; b >>= 1, a = 1ll * a * a % mod)
        if (b & 1) c = 1ll * c * a % mod;
    return c % mod;
}


/**
 * Fast Power
 * @return a^b
 * @timecomplex O(nlgn)
 */
MATHTYPE fp(MATHTYPE a, MATHTYPE b) {
    MATHTYPE c = 1;
    for (; b; b >>= 1, a = 1ll * a * a)
        if (b & 1) c = 1ll * c * a;
    return c;
}


/**
 * Fast Multiply with mod
 * @param a,b < 2^63-1
 * @return a*b%m
 * @timecomplex O(algb)
 */
MATHTYPE multiMod(MATHTYPE a, MATHTYPE b, MATHTYPE mod) {
    b %= mod;
    a %= mod;
    MATHTYPE res = 0;
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


/**
 * Judge whether a*b overflows
 */
//判断两数相乘是否越界
MATHTYPE valid_multi(MATHTYPE a, MATHTYPE b) {
    MATHTYPE MAX = (MATHTYPE)1e18, ret = 0;
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


/**
 * atoi with mod
 */
MATHTYPE atoi_mod(string val, MATHTYPE mod) {
    MATHTYPE res = 0;
    for (MATHTYPE i = 0; i < val.length(); ++i) {
        res = ((res)*10 + val[i] - '0') % mod;
    }
    return res;
}


/**
 * Euclid Algorithm
 * @return great common divider of a and b
 */
MATHTYPE gcd(MATHTYPE a, MATHTYPE b) {
    if (a < b) swap(a, b);
    while (b != 0) {
        MATHTYPE c = a;
        a = b;
        b = c % b;
    }
    return a;
}


/**
 * LCM
 */
MATHTYPE lcm(MATHTYPE a, MATHTYPE b) {
    MATHTYPE g = gcd(a, b);
    return a / g * b;
}


/**
 * Extended Euclid Algorithm
 * @return gcd of a and b
 * @return x,y a*x+b*y=gcd(a,b)
 */
MATHTYPE ext_gcd(MATHTYPE a, MATHTYPE b, MATHTYPE& x, MATHTYPE& y) {
    if (b == 0) {
        x = 1;
        y = 0;
        return a;
    } else {
        MATHTYPE r = ext_gcd(b, a % b, y, x);
        y -= x * (a / b);
        return r;
    }
}


/**
 * Inversion of mod
 * The first one uses Extended Euclid Algorithm
 * (Deprecated)The second one uses Fast Mod. ONLY FOR prime.
 */
MATHTYPE inv(MATHTYPE a, MATHTYPE mod) {
    a %= mod;
    MATHTYPE x, y;
    ext_gcd(mod, a, x, y);
    return md(y);
}
MATHTYPE inv2(MATHTYPE a, MATHTYPE mod){
    return fp(a,mod-2,mod);
}


/**
 * Prepare the factorial and its inversion
 * @timecomplex O(n)
 */
MATHTYPE fac[MAXN], fac_inv[MAXN];
void pre_fac_inv(MATHTYPE n) {
    fac[0] = 1;
    for (MATHTYPE i = 1; i <= n; ++i) {
        fac[i] = md(1ll * fac[i - 1] * i);
    }
    fac_inv[n] = inv(fac[n], mod);
    for (MATHTYPE i = n - 1; i >= 0; ++i) {
        fac[i] = md(1ll * fac[i + 1] * (i + 1));
    }
}


/**
 * Prepare the first n inversion and factorial inversion
 * @timecomplex O(n)
 */
MATHTYPE n_inv[MAXN];
void pre_n_inv(MATHTYPE n) {
    n_inv[1] = 1;
    for (MATHTYPE i = 2; i <= n; i++) {
        n_inv[i] = md(1ll * md(-mod / i) * md(n_inv[mod % i]));
    }

    fac_inv[0] = fac_inv[1] = 1;
    for (MATHTYPE i= 2;i<=n;i++){
        fac_inv[i] = fac_inv[i-1]*n_inv[i];
    }
}


/**
 * Lucas Agorithm
 * @param mod mod should be prime
 * @return C(n,m)
 * @timecomplex O(nlgp)
 */
MATHTYPE lucas(MATHTYPE n, MATHTYPE m, MATHTYPE mod) {
    if (n < m)
        return 0;
    else if (n < mod)
        return md(md(fac[n] * fac_inv[m]) * fac_inv[n - m]);
    else
        return md(lucas(n / mod, m / mod, mod) * lucas(n % mod, m % mod, mod));
}


/**
 * Linear Time Prepare Prime
 * @timecomplex O(n)
 */
MATHTYPE prime_table(MATHTYPE n) {
    bool vis[MAXN];
    MATHTYPE prime[MAXN];
    MATHTYPE used = 0;
    vis[0] = vis[1] = 1;
    for (MATHTYPE i = 2; i <= n; ++i) {
        if (!vis[i]) {
            prime[used++] = i;
            //积性函数f(p)
            // phi(p) = p-1
        }
        for (MATHTYPE j = 0; j < used; ++j) {
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


/**
 * Segment Prime Preparation
 * Based on Egypt Sieve Algorithm
 * @timecomplex O(nlgnlgn)
 */
void segment_prime(MATHTYPE a, MATHTYPE b, bool p[], bool sqrtbp[]) {
    for (MATHTYPE i = 2; i * i <= b; ++i)
        if (!sqrtbp[i]) {
            for (MATHTYPE j = i * i; j * j <= b; j += i) sqrtbp[j] = 1;
            for (MATHTYPE j = max(i * i, (a - 1) / i + 1) * i; j <= b; j += i)
                p[j - a] = 1;
        }
}


/**
 * Euler Function
 * @return phi(n)
 */
MATHTYPE phi(MATHTYPE n) {
    MATHTYPE ans = n;
    for (MATHTYPE i = 2; i * i <= n; ++i)
        if (ans % i == 0) {
            ans = ans / i * (i - 1);
            while (n % i == 0) n /= i;
        }
    if (n > 1) ans = ans / n * (n - 1);
    return ans;
}


/**
 * Euler Power Reduction Formula
 * @return a^b%m
 */
MATHTYPE FM_Euler(MATHTYPE a, MATHTYPE b, MATHTYPE m) {
    MATHTYPE p = phi(m);
    b = b % p;
    return fp(a, b, m);
}


/**
 * Solve Congruence Equation
 * @return ans for equation a mod r = b
 */
bool cong_eq(MATHTYPE a, MATHTYPE b, MATHTYPE m, MATHTYPE ans[], MATHTYPE& len) {
    MATHTYPE d, x, y;
    d = ext_gcd(a, m, x, y);
    if (b % d) return false;
    MATHTYPE base = ((b / d * x) % m + m) % m;
    len = d;
    for (MATHTYPE i = 0; i < len; ++i) ans[i] = (base + i * (m / len)) % m;
    return true;
}


/**
 * China Remainer Theorem
 * Solve Congruence Equation Set
 * @param b[] mod should be prime
 * @return ans of equation set a mod r = b
 */
MATHTYPE china(MATHTYPE a[], MATHTYPE b[], MATHTYPE r) {
    // x = a[i] mod b[i];
    MATHTYPE M = 1, i, Mi, ans = 0;
    for (i = 0; i < r; ++i) M *= b[i];
    for (i = 0; i < r; ++i) {
        Mi = M / b[i];
        ans = (ans + inv(Mi, b[i]) * Mi * a[i]) % M;
    }
    return (ans + M) % M;
}


/**
 * Solve Congruence Equation Set
 */
MATHTYPE cong_eq_group(MATHTYPE a[], MATHTYPE b[], MATHTYPE r) {
    MATHTYPE ans = 0, M = 1, x, y, d;
    for (MATHTYPE i = 0; i < r; ++i) {
        d = ext_gcd(M, b[i], x, y);  // b[i]*y=d(mod M)
        if ((ans - a[i]) % d) return -1;
        M = M / d * b[i];
        MATHTYPE z = y * ((ans - a[i]) / d);
        ans = z * b[i] + a[i];
        ans = ((ans % M) + M) % M;
    }
    return ans;
}


/**
 * Discrete logarithm
 * @return t st. A^t % mod = C
 */
MATHTYPE BSGSx(MATHTYPE A, MATHTYPE C, MATHTYPE mod) {
    A %= mod;
    C %= mod;
    if (C == 1) return 0;
    MATHTYPE cnt = 0, tmp = 1;

    for (MATHTYPE g = gcd(A, mod); g != 1; g = gcd(A, mod)) {
        if (C % g) return -1;
        C /= g;
        mod /= g;
        tmp = tmp * A / g % mod;
        ++cnt;
        if (C == tmp) return cnt;
    }
    // tmp*a^(x-cnt)=b' (mod c')

    MATHTYPE T = (MATHTYPE)sqrt(0.5 + mod);
    MATHTYPE b = C;
    std::map<MATHTYPE, MATHTYPE> hash;
    hash[b] = 0;
    for (MATHTYPE i = 1; i <= T; ++i) {
        b = b * A % mod;
        hash[b] = i;
    }

    A = fp(A, T, mod);
    for (MATHTYPE u = 1; u <= T; ++u) {
        tmp = tmp * A % mod;
        if (hash.count(tmp)) return u * T - hash[tmp] + cnt;
    }
    return -1;
}


/**
 * MillerRabin Prime Test
 * @param times test times, the probability of mistake is 4^(-times)
 * @return whether n is a prime
 */
bool witness(MATHTYPE s, MATHTYPE n) {
    MATHTYPE u = n - 1;
    MATHTYPE t = 0;
    while ((u & 1) == 0) u >>= 1, t++;  //Second detection theorem + Fermat Test

    MATHTYPE x = fp(s, u, n), tmp;
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

bool millerRabin(MATHTYPE n, const MATHTYPE times = 3) {
    if (n == 2) return true;
    if ((n & 1) == 0 || n < 2) return false;
    MATHTYPE i = times;
    while (i--) {
        MATHTYPE s = randint(1, n - 1);
        if (witness(s, n)) return false;  // s^(n-1)==1 (mod n)
    }
    return true;
}


/** 
 * Pollard Decomposition of prime factors
 * @timecomplex O(n^(1/4))
 */
MATHTYPE pollard_rho(MATHTYPE n) {
    MATHTYPE x, y, k = 2, d, i = 1, c;
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

void split(MATHTYPE n, MATHTYPE factors[], MATHTYPE& len) {
    if (millerRabin(n))
        factors[++len] = n;  // or ++fac[n] or (map<MATHTYPE,MATHTYPE> ++fac[n])
    else {
        MATHTYPE p;
        do {
            p = pollard_rho(n);
        } while (p == n || p == 1);
        split(p, factors, len);
        split(n / p, factors, len);
    }
}


/**
 * Decomposition of prime factors
 * @timecomplex O(sqrt(n)/logn)
 * @return number of prime factors
 */
MATHTYPE prime[MAXN];
MATHTYPE split_(MATHTYPE n, MATHTYPE a[], MATHTYPE& len) {
    MATHTYPE sqrtn = (MATHTYPE)sqrt(n);
    len = 0;
    a[0] = 1;
    MATHTYPE ans = 1, tmp;
    for (MATHTYPE i = 0; prime[i] <= sqrtn; ++i) {
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


/**
 * Order of a number
 */
MATHTYPE ord(MATHTYPE a, MATHTYPE m) {
    MATHTYPE p = phi(m), tmp = 1;
    for (MATHTYPE i = 1; i <= p; i++) {
        tmp = tmp * a % m;
        if (tmp == 1) return i;
    }
    return -1;
}


/**
 * Primitive root
 */
MATHTYPE is_g(MATHTYPE x, MATHTYPE factors[], MATHTYPE& len, MATHTYPE phi) {
    for (MATHTYPE i = 0; i < len; ++i)
        if (factors[i] != 1 && fp(x, phi / factors[i], phi) == 1) return false;
    return true;
}

MATHTYPE get_g(MATHTYPE p) {
    if (p == 2 || p == 4) return 1;
    if (!(millerRabin(p) || (!(p & 1) && millerRabin(p / 2)))) return -1;
    MATHTYPE phiofn = phi(p), factors[1000], len = 0;
    split(phiofn, factors, len);
    for (MATHTYPE i = 1; i < p; i++)
        if (is_g(i, factors, len, phiofn)) return i;
    return -1;
}


/**
 * Number of Relatively Prime Number with a in [1,n]
 * -First way uses inclusion and exclusion principle
 * -Second way uses Mobius inversion
 */
MATHTYPE prime_of_n[100];
MATHTYPE mutual_prime_(MATHTYPE a, MATHTYPE n) {
    MATHTYPE len;
    split(a, prime_of_n, len);
    MATHTYPE mask = 1 << len, ans = 0;
    for (MATHTYPE i = 1; i < mask; ++i) {
        MATHTYPE cnt = 0, v = 1;
        for (MATHTYPE j = 0; j < len; ++j)
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

MATHTYPE mu[100];
MATHTYPE mutual_prime(MATHTYPE a, MATHTYPE b) {
    MATHTYPE ans;
    if (a > b) std::swap(a, b);
    for (MATHTYPE i = 1; i <= a; ++i) ans += mu[i] * (a / i) * (b / i);
    return ans;
}

#endif

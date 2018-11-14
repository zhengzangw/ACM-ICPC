#ifndef STRUCTRUE_H
#define STRUCTRUE_H

/**
 * Head Info
 */
#include <bits/stdc++.h>
using namespace std;
typedef long long DTYPE;
using namespace std;
#define MAXN 0
int n;
DTYPE a[MAXN]; //@index 1

/**
 * Prefix Sum
 * @index 1
 */
DTYPE presum[MAXN];
void presum_init(int n, DTYPE* a){
    presum[1] = a[1];
    for (int i=2;i<=n;++i){
        presum[i] = presum[i-1] + a[i];
    }
}


/**
 * Discrete
 * @param n number of elements
 * @param a array to be discreted
 * @return len number of different elements in a
 * @return val_a val_a[i] means the ith smallest elememt in a
 * @return a dis_a[a[i]] is the original a
 */
DTYPE dis_a[MAXN], val_a[MAXN];
DTYPE discrete(int n, const DTYPE a[]){
    for (int i=0;i<n;++i){
        dis_a[i] = val_a[i] = a[i];
    }
    sort(val_a,val_a+n);
    int len = int(unique(val_a,val_a+n) - val_a);
    for (int i=0;i<n;++i) dis_a[i] = int(lower_bound(val_a, val_a+len, a[i])- val_a + 1);
    return len;
}


/***
 * Interval Function
 * Interval Union Function of O(f(n)) F(n) = G(F(k),F(n-k))
 * Interval Add Function of O(f(n)) F(n) = G(F(k)+F(n-k))
 * Interval Difference Function of O(f(n)) F(n) = G(F(k),F(k-n))
 * Interval Minus Funtion of O(f(n)) F(n) = G(F(k)-F(k-n))
 */

/**
 * Binary Indexed Tree
 * @op Modify(Singal) O(lgn)
 * @op Query(Interval, MinusF) O(lgn)
 * ind x contains from lowbit(x) to x
 * @index start from 1
 */
DTYPE BIT[MAXN];
inline static DTYPE lowbit(DTYPE x){
    return x&-x;
}

inline void BIT_add(int ind, DTYPE k){
    while (ind<=n){
        BIT[ind] += k;
        ind += lowbit(ind);
    }
}

inline DTYPE BIT_presum(int ind){
    DTYPE ret = 0;
    while (ind){
        ret += BIT[ind];
        ind -= lowbit(ind);
    }
    return ret;
}

inline void BIT_update(int ind, DTYPE k){
    while (ind<=n){
        BIT[ind] = min(BIT[ind],k);
        ind += lowbit(ind);
    }
}

inline DTYPE BIT_premin(int ind){
    DTYPE min_val = INT_MAX;
    while (ind){
        min_val = min(min_val, BIT[ind]);
        ind -= lowbit(ind);
    }
    return min_val;
}

// With differential, can realize interval update but single point lookup
inline void BIT_init(int n, const DTYPE *a){
    memset(BIT,0,sizeof(DTYPE)*(n+1));
    for (int i=1;i<=n;++i){
        BIT[i] += a[i];
        if (i+lowbit(i)<=n)
            BIT[i+lowbit(i)] += BIT[i];
    }
}



#endif
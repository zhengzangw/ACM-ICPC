#include <iostream>
#include <cstdio>
#include <queue>
using namespace std;
typedef unsigned long long ull;
const ull Base=100000007; //1844674407370955161
//const ull prime = 1844674407370955161; 自然溢出
// 由于冲突概率低，忽略相等时再次比较
// 生日攻击理论： 对于哈希值在[0,n)均匀分布的哈希函数，首次冲突发生的期望是O(sqrt(n))
int n,m,t,cnt;


int contain(string b, string a)
{
    int a1 = a.length(), b1 = b.length();
    if (a1>b1) return false;
    
    ull t = 1;
    for (int i=0; i<a1; i++) t*=Base;
    
    ull ah = 0, bh = 0;
    for (int i=0;i<a1;i++) ah = ah*Base + a[i];
    for (int i=0;i<a1;i++) bh = bh*Base + b[i];
    for (int i = 0; i+a1 <= b1; i++){
        if (ah == bh) return i;
        if (i+a1<b1) bh = bh * Base + b[i + a1] - b[i]*t;
    }
    return -1;
}

int overlap(string b, string a)
{
    int a1 = a.length(), b1 = b.length();
    int ans = 0;
    ull ah = 0, bh = 0, t =1;
    for (int i = 1; i<= min(a1,b1); i++)
    {
        ah = ah + a[a1-i]*t;
        bh = bh*Base + b[i-1];
        if (ah == bh) ans = i;
        t *= B;
    }
    return ans;
}

int main(int argc, char *argv[])
{
    string s,t;
    cin >> s >> t;
    cout << contain(s,t)+1;
    
    return 0;
}

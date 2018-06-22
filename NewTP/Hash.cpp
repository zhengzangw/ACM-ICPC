#include <iostream>
#include <cmath>
#include <cstdio>
#include <set>
#include <map>
#include <unordered_map>
#define fabs(x) ((x) > 0 ? (x) : -(x))
using namespace std;

class HashOA
{
  private:
    struct node
    {
        int key, cnt;
    };
    static const int m = 8000037; //Memory 8*m/10^6 MB Time O(n/(m-n))
    static const int mod = 4000037;
    node Hash[m];
    inline int f_division(int k)
    {
        if (k < 0)
            return (k % m) + m;
        return k % m;
    }
    inline int f_division2(int k)
    {
        if (k < 0)
            return (k % mod) + mod;
        return k % mod;
    }
    int f_multi(int k)
    {
        static const double A = 0.6180339887498949;
        double t = abs(k * A - ceil(k * A));
        return ceil(m * t);
    }
    inline int linearprobe(int k, int i)
    {
        return (f_division(k) + i) % m;
    }
    inline int quadraticprobe(int k, int i)
    {
        return (f_division(k) + ((i + 1) * i) / 2) % m;
    }
    inline int doubleprobe(int k, int i)
    {
        return (f_division(k) + i * f_division2(k)) % m;
    }

  public:
    int insert(int k)
    {
        int i=0,j = f_division(k);
        while (i++ < m)
        {
            if (Hash[j].cnt <= 0 || Hash[j].key == k)
            {
                Hash[j].key = k;
                Hash[j].cnt++;
                return j;
            }
            j = doubleprobe(k,i);
        }
        return -1;
    }
    int search(int k)
    {
        int i=0,j = f_division(k);
        while (i++ < m && Hash[j].cnt != 0)
        {
            if (Hash[j].key == k)
                return Hash[j].cnt;
            j = doubleprobe(k, i);
        }
        return 0;
    }
    int del(int k)
    {
        int i = 0, j = f_division(k);
        while (i++ < m && Hash[j].cnt != 0)
        {
            if (Hash[j].key == k){
                Hash[j].cnt --;
                return j;
            }
            j = doubleprobe(k, i);
        }
        return -1;
    }
};
HashOA H;


int main()
{
    //Use set as Hash table
    //Based on BST
    set<int> s;
    s.insert(1);

    //Use map as Hash table (can store occur time)
    //Based on BST O(lgn)
    map<int,int> h;
    h[1] = 1;
    //中序遍历BST即可获得排好序的key列表
    //BST可以很容易地执行顺序统计，找出最接近的最小和最大元素

    //Based on Hash table
    unordered_map<int,int> Hash;
    Hash[1] = 1;

    /********Some hashing method**********
    >Radix presentation
    >Sum
    >Sum of spuare
    >xor
    >v[i]<<(i*3) xor
    >Multiplication*Prime
    >(PartialSum)^|&(PartialSum)
    */
    return 0;
}

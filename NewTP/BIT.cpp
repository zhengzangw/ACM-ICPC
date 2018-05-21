#include <iostream>
#include <vector>
using namespace std;
class BIT : vector<int>
{
  private:
    inline int lowbit(int k) { return k & -k; }

  public:
    explicit BIT(int k = 0)
    {
        assign(k + 1, 0);
    }
    int getsum(int k)
    {
        return k > 0 ? getsum(k - lowbit(k)) + (*this)[k] : 0;
    }
    int getsum(int L, int R)
    {
        return getsum(R) - getsum(L - 1);
    }
    int last() { return int(size()) - 1; };
    void update(int k, int u)
    {
        while (k <= last())
        {
            (*this)[k] += u;
            k += lowbit(k);
        }
    }
};
class BITadvanced : vector<int>
{
  private:
    inline int lowbit(int k) { return k & -k; }
    vector<int> lazy;
    void update_lazy(int k, int u)
    {
        while (k <= last())
        {
            lazy[k] += u;
            k += lowbit(k);
        }
    }
    int getsum(int k, vector &a)
    {
        return k > 0 ? getsum(k - lowbit(k), a) + a[k] : 0;
    }

  public:
    explicit BITadvanced(int k = 0)
    {
        assign(k + 1, 0);
        lazy.assign(k + 1, 0);
    }

    int getsum(int L, int R)
    {
        int s = 0;
        s += getsum(R, *this) + getsum(R, lazy) * R;
        s -= getsum(L - 1, *this) + getsum(L - 1, lazy) * (L - 1);
        return s;
    }
    int last() { return int(size()) - 1; };
    void update(int k, int u)
    {
        while (k <= last())
        {
            (*this)[k] += u;
            k += lowbit(k);
        }
    }
    void update(int L, int R, int u)
    {
        update(L, -(L - 1) * u);
        update_lazy(L, u);
        update(R + 1, R * u);
        update_lazy(R + 1, -u);
    }
};

class BIToptimal : vector<int> //low efficiency
{
  private:
    inline int lowbit(int k) { return k & -k; }
    int *num;

  public:
    explicit BIToptimal(int k = 0, int *a = NULL)
    {
        num = a;
        assign(k + 1, -INT_MAX);
        for (int i = 1; i <= k; i++)
        {
            (*this)[i] = num[i];
            for (int j = 1; j < lowbit(i); j <<= 1)
                (*this)[i] = max((*this)[i], (*this)[i - j]);
        }
    }
    int getmax(int L, int R)
    {
        int ans = num[R];
        while (L != R)
        {
            for (R -= 1; R - L >= lowbit(R); R -= lowbit(R))
                ans = max(ans, (*this)[R]);
            ans = max(ans, num[R]);
        }
        return ans;
    }
    int last() { return int(size()) - 1; };
    void update(int k, int u)
    {
        num[k] += u;
        for (int i = k; i <= last(); i += lowbit(i))
        {
            (*this)[i] = u;
            for (int j = 1; j < lowbit(i); j <<= 1)
            {
                (*this)[i] = max((*this)[i], (*this)[i - j]);
            }
        }
    }
};

int n, m, t, x, y;
int a[51000], b[51000];
int main()
{
    cin >> n >> m;
    BITadvanced bit(n);
    for (int i = 1; i <= n; i++)
    {
        cin >> t;
        bit.update(i, t);
    }
    for (int i = 1; i <= m; i++)
    {
        cin >> x >> y;
        bit.update(x, y, 1);
        cin >> x >> y;
        cout << bit.getsum(x, y) << endl;
    }
    return 0;
}

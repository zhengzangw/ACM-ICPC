#include <iostream>
#include <vector>
using namespace std;
const int MAX_N = 1 << 20;
int par[MAX_N], ran[MAX_N];
int n, m, t;
char c;

class DisjointSet : vector<int>
{
  private:
    vector<int> rank;
    vector<int> num;
    int find(int x)
    {
        if ((*this)[x] == x)
            return x;
        else
            return ((*this)[x]) = find((*this)[x]);
    }

  public:
    explicit DisjointSet(int k = 0)
    {
        assign(k + 1, 0);
        for (int i = 1; i <= k; i++)
            (*this)[i] = i;
        rank.assign(k + 1, 0);
        num.assign(k + 1, 1);
    }
    int same(int x, int y) { return find(x) == find(y); }
    void merge(int x, int y)
    {
        x = find(x);
        y = find(y);
        if (x == y)
            return;
        if (rank[x] < rank[y])
        {
            (*this)[x] = y;
            num[y] += num[x];
        }
        else
        {
            (*this)[y] = x;
            num[x] += num[y];
        }
        if (rank[x] == rank[y])
            rank[x]++;
    }
    int getsetlen(int x) { return num[find(x)]; }
};

int main()
{
    cin >> n;
    DisjointSet s(n);
    for (int i = 0; i < n; i++)
    {
        cin >> c >> m >> t;
        if (c == 'q')
        {
            if (s.same(m, t))
                cout << "Yes\n";
            else
                cout << "No\n";
        }
        if (c == 'u')
        {
            s.merge(m, t);
            cout << "Done!\n";
        }
        if (c == 'n')
        {
            cout << s.getsetlen(m) << ' ' << s.getsetlen(t) << endl;
        }
    }
}

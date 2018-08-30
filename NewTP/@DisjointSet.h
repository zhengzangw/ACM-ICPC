#ifndef DisjointSet_H
#define DisjointSet_H
#include <vector>

class DisjointSet : std::vector<int>
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

#endif

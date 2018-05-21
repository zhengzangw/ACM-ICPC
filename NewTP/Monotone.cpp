#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iomanip>
typedef long long LL;
using namespace std;
#define MAX_N 1000
struct Node
{
    int key, limit;
    Node(int _key = 0, int _limit = 0)
    {
        key = _key;
        limit = _limit;
    }
    bool operator>(const Node &a) const
    {
        return (key > a.key);
    }
};
class MonotoneQueue //不严格单调增
{
  private:
    int limit, head, tail;
    Node Q[MAX_N];

  public:
    MonotoneQueue()
    {
        head = tail = 0;
    }
    bool push(Node x)
    {
        if (head < tail && x.limit < Q[tail].limit)
            return false;
        while (head < tail && Q[tail - 1] > x)
            tail--;
        Q[tail++] = x;
        return true;
    }
    bool push(int key, int limit)
    {
        return push(Node(key, limit));
    }
    void updatelimit(int _limit) { limit = _limit; }
    bool top(Node &rt)
    {
        while (tail > head && Q[head].limit < limit)
            head++;
        if (head == tail)
            return false;
        rt = Q[head];
        return true;
    }
};

class MonotoneStack
{
  private:
    int Q[MAX_N];
    int tail;

  public:
    MonotoneStack()
    {
        tail = 0;
    }
    void push(int x)
    {
        while (tail > 0 && Q[tail - 1] > x)
            tail--;
        Q[tail++] = x;
    }
    bool top(int &rt)
    {
        if (tail == 0)
            return false;
        else
            rt = Q[tail];
        return true;
    }
    bool pop()
    {
        if (tail == 0)
            return false;
        else
            tail--;
        return true;
    }
    void print()
    {
        for (int i = 0; i < tail; i++)
        {
            cout << Q[i];
        }
        cout << endl;
    }
};

int main()
{
    MonotoneStack S;
    int n;
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        cin >> i;
        S.push(i);
        S.print();
    }
    return 0;
}

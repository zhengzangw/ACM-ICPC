#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <vector>
typedef long long LL;
using namespace std;
class Haffman
{
  private:
    struct Node
    {
        int weight, flag, parent, lch, rch;
    };
    struct Code
    {
        int bit[6];
        int len;
        int weight;
    };
    static const int MAX_N = 100010;
    int n, TotalWPL;
    Node HT[MAX_N];
    Code HTcode[MAX_N];
    void combine(int x1, int x2, int i)
    {
        HT[x1].parent = n + i;
        HT[x2].parent = n + i;
        HT[x1].flag = 1;
        HT[x2].flag = 1;
        HT[n + i].weight = HT[x1].weight + HT[x2].weight;
        HT[n + i].lch = x1;
        HT[n + i].rch = x2;
        TotalWPL += HT[x1].weight + HT[x2].weight;
    }

  public:
    Haffman(int len, int w[])
    {
        n = len;
        TotalWPL = 0;
        int m1 = 0, m2 = 0, x1 = 0, x2 = 0;
        for (int i = 0; i < 2 * n - 1; i++)
        {
            if (i < n)
            {
                HT[i].weight = w[i];
            }
            else
                HT[i].weight = 0;
            HT[i].flag = 0;
            HT[i].parent = HT[i].lch = HT[i].rch = -1;
        }
        for (int i = 0; i < n - 1; i++)
        {
            m1 = m2 = INT_MAX;
            for (int j = 0; j < n + i; j++)
            {
                if (!HT[j].flag && HT[j].weight < m1)
                {
                    m2 = m1;
                    x2 = x1;
                    m1 = HT[j].weight;
                    x1 = j;
                }
                else if (!HT[j].flag && HT[j].weight < m2)
                {
                    m2 = HT[j].weight;
                    x2 = j;
                }
            }
            combine(x1, x2, i);
        }
    }
    int WPL(int x)
    {
        int sum = 0;
        for (Node i = HT[x]; i.parent != -1; i = HT[i.parent])
            sum++;
        return sum * HT[x].weight;
    }
    int WPL()
    {
        return TotalWPL;
    }
    void generateHFcode()
    {
        Code *cd = new Code;
        int child, parent;
        for (int i = 0; i < n; i++)
        {
            cd->len = 0;
            cd->weight = HT[i].weight;
            child = i;
            parent = HT[child].parent;
            while (parent != -1)
            {
                if (HT[parent].lch == child)
                    cd->bit[cd->len] = 0;
                else
                    cd->bit[cd->len] = 1;
                cd->len++;
                child = parent;
                parent = HT[child].parent;
            }
            for (int j = cd->len; j > 0; j--)
                HTcode[i].bit[cd->len - j] = cd->bit[j - 1];
            HTcode[i].len = cd->len;
            HTcode[i].weight = cd->weight;
        }
    }
    int getcodeweight(int i)
    {
        return HTcode[i].weight;
    }
    string getcode(int i)
    {
        char b[6] = "\0";
        for (int j = 0; j < HTcode[i].len; j++)
            b[j] = HTcode[i].bit[j] + '0';
        b[HTcode[i].len] = '\0';
        string s = b;
        return s;
    }
    string encode(string s, int dict[])
    {
        string ans;
        for (int i = 0; i < s.length(); i++)
        {
            ans += getcode(dict[s[i]]);
        }
        return ans;
    }
    string decode(string s, char revdict[])
    {
        string ans;
        int node;
        for (int i = 0; i < s.length();)
        {
            node = 2 * n - 2;
            while (HT[node].lch != -1)
            {
                if (s[i] == '0')
                    node = HT[node].lch;
                else
                    node = HT[node].rch;
                i++;
            }
            ans += revdict[node];
        }
        return ans;
    }
};

int main()
{
    int a[100], n, revdict[1000];
    char dict[1000];
    char c;
    string s;
    cin >> n;
    for (int i = 0; i < n; i++)
        cin >> a[i];
    for (int i = 0; i < n; i++)
    {
        cin >> c;
        revdict[c] = i;
        dict[i] = c;
    }
    cin >> s;
    Haffman h(n, a);
    cout << h.WPL() << endl;
    h.generateHFcode();
    for (int i = 0; i < n; i++)
    {
        cout << dict[i] << ": Weight=" << h.getcodeweight(i) << " Code=" << h.getcode(i) << endl;
    }
    cout << h.encode(s, revdict) << endl;
    cout << h.decode(h.encode(s, revdict), dict) << endl;
    return 0;
}

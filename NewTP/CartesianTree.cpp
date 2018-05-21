#include <iostream>
#include <cstdio>
#include <string>
#include <cstring>
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <stack>
typedef long long LL;
using namespace std;
class CartesianTree
{
  private:
	static const int MAX_N = 100010;
	struct node
	{
		int par, lch, rch;
		LL val;
	};
	node T[MAX_N];
	int N, root;
	LL ans;
	int dfs(int cur)
	{
		if (cur == -1)
			return 0;
		int num = dfs(T[cur].lch) + dfs(T[cur].rch) + 1;
		ans = max(ans, num * T[cur].val);
		return num;
	}

  public:
	//小根
	CartesianTree(LL a[], int n)
	{
		N = n;
		int st[MAX_N], k, top = -1;
		// 右链更新法
		for (int i = 0; i < n; i++)
		{
			T[i].val = a[i];
			T[i].lch = -1;
			T[i].rch = -1;
			T[i].par = -1;
			k = top;
			while (k >= 0 && a[st[k]] > a[i])
				k--;
			if (k != -1)
			{
				T[i].par = st[k];
				T[st[k]].rch = i;
			}
			if (k < top)
			{
				T[st[k + 1]].par = i;
				T[i].lch = st[k + 1];
			}
			st[++k] = i;
			top = k;
		}
		root = st[0];
	}
	LL DFS()
	{
		ans = 0;
		dfs(root);
		return ans;
	}
};
int main()
{
	LL a[100010];
	int n;
	while (cin >> n && n != 0)
	{
		for (int i = 0; i < n; i++)
			cin >> a[i];
		CartesianTree T(a, n);
		cout << T.DFS() << endl;
	}
	return 0;
}

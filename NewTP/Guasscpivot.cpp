#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <stack>
#include <iomanip>
#include <cstdlib>
typedef long long LL;
typedef unsigned int UI;
typedef unsigned long long ULL;
const LL mod = 1000000007;
#define MAXN 100
#define eps 1e-10
#define fabs(x) ((x) > 0 ? (x) : -(x))
using namespace std;

int gauss_cpivot(int n, double a[][MAXN], double b[])
{
    // Solve A*X=b, Ans is in b[] O(n^3)
    int i, j, k, row;
    double maxp, t;
    // For every column
    for (k = 0; k < n; k++)
    {
        // Find the max element
        for (maxp = 0, i = k; i < n; i++)
            if (fabs(a[i][k]) > fabs(maxp))
                maxp = a[row = i][k];
        // If failed, then there is no only solution
        if (fabs(maxp) < eps)
            return 0;
        // Exchange Row row and Row k
        if (row != k)
        {
            for (j = k; j < n; j++)
                t = a[k][j], a[k][j] = a[row][j], a[row][j] = t;
            t = b[k], b[k] = b[row], b[row] = t;
        }
        // Update rest elements
        for (j = k + 1; j < n; j++)
        {
            a[k][j] /= maxp;
            for (i = k + 1; i < n; i++)
                a[i][j] -= a[i][k] * a[k][j];
        }
        b[k] /= maxp;
        for (i = k + 1; i < n; i++)
            b[i] -= b[k] * a[i][k];
    }
    // Deal with the 上三角矩阵
    for (i = n - 1; i >= 0; i--)
        for (j = i + 1; j < n; j++)
            b[i] -= a[i][j] * b[j];
    return 1;
}

int main()
{
    //std::ios::sync_with_stdio(false);
    int m, n, t, ans;
    double a[MAXN][MAXN], b[MAXN];
    cin >> n;
    for (int i = 0; i < n; i++)
    {
        for (int j = 0; j < n; j++)
            cin >> a[i][j];
        cin >> b[i];
    }
    if (gauss_cpivot(n, a, b))
        for (int i = 0; i < n; i++)
            cout << "X" << i + 1 << "=" << b[i] << endl;
    else
        cout << "No solution or No only solution!";

    return 0;
}

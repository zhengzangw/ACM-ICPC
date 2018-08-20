#ifndef header_H
#define header_H

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
#define mod 1000000007;
#define eps 1e-10
#define max(a, b) (((a) > (b)) ? (a) : (b))
#define min(a, b) ((a)>(b))? (b) : (a))
#define fabs(x) ((x) > 0 ? (x) : -(x))
#define randint(x, y) rand() % (y) + (x)
using int64 = long long;
using namespace std;

inline int readint() //快读
{
    int ret(0);
    int sgn(1);
    char c;
    while (isspace(c = getchar()));
    if (c == '-')
    {
        sgn = -1;
        c = getchar();
    }
    while (isdigit(c))
    {
        ret = (ret << 3) + (ret << 1) + c - '0';
        c = getchar();
    }
    return ret * sgn;
}

#endif
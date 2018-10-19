#ifndef Graph_H
#define Graph_H

#define MAXN 0
#include <algorithm>
using namespace std;

/**
 * Edge
 * @type Data Structure
 */
struct Edge{
    int u,v,w;
    bool operator<(const Edge &t) const
    {
        return w<t.w;
    }
};
Edge edge_set[3*MAXN]; // start index from 1

/**
 * Disjoint Set
 * @type Data Structure
 */
int rnk[MAXN]; //Merge with rank
int DJS[MAXN]; //Disjoint set
int num[MAXN]; //Attr: num, max, min, rel(or use factDJS)

static int find(int x)
{
    if (DJS[x]==x) return x;
    else return (DJS[x]=find(DJS[x]));
}

int eq(int x, int y)
{
    return (find(x)==find(y));
}

// if weighted-DJS is used, do not use rank optimization
void merge(int x, int y)
{
    x = find(x);
    y = find(y);
    if (x==y) return;
    if (rnk[x] < rnk[y]){
        DJS[x] = y;
        num[y] = num[x] + num[y]; //Attr alters
    } else {
        DJS[y] = x;
        num[x] = num[y] + num[x]; //Attr alters
    }
    if (rnk[x] == rnk[y])
        rnk[x]++;
}

void DJS_initial(){
    for (int i=0;i<MAXN;++i){
        DJS[i] = num[i] = i;
        rnk[i] = 0;
    }
}

/**
 * Kruskal 
 * @type Algotithm for minimal spinning tree
 * @param edge_set
 * @param n vertices number
 * @param m edge number
 * @return weight of minimal spinning tree
 * @timecomplex O(Elog(V))
 */
int Kruskal(int n, int m){
    int ret = 0;
    DJS_initial();
    sort(edge_set+1, edge_set+n+1);
    for (int i=1;i<m;++i)
        if (!eq(edge_set[i].u,edge_set[i].v)) {
            merge(edge_set[i].u,edge_set[i].v);
            ret += edge_set[i].w;         
        }
    // Lack: Judge whether the graph is connected
    return ret;
}

#endif

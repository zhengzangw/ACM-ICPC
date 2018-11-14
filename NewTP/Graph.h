#ifndef Graph_H
#define Graph_H

#include <bits/stdc++.h>
#include "math.h"
#include "matrix.h"
using namespace std;
#define MAXV 0
#define MAXE 0
typedef int WTYPE;
typedef pair<WTYPE, int> pTi;

/**
 * Edge
 * @type Data Structure
 */
struct Edge {
    int u, v;
    WTYPE w;
    bool operator<(const Edge &t) const { return w < t.w; }
};
Edge edge_set[MAXE];  // start index from 1
int edge_head;

/**
 * Adjacency Matrix
 * @type Data Structure
 */
WTYPE adj_mtrx[MAXV][MAXV];
int deg[MAXV];

/**
 * Adjacency List
 * @type Data Structure
 */
struct Edgenode {
    int v;
    WTYPE w;
    Edgenode *nxt;
};
struct Vertexnode {
    int etotal;
    Edgenode *head;
    Vertexnode() : etotal(0), head(NULL) {}
};
struct Adjlist {
    Vertexnode v[MAXV];
    int vtotal, etotal;
    Adjlist() : vtotal(0), etotal(0) {}
} adj_list;
// add u->v
void add_edge_adjlist(int u, int v, WTYPE w) {
    Edgenode *p = new Edgenode;
    p->v = v;
    p->w = w;
    p->nxt = adj_list.v[v].head;
    adj_list.v[v].head = p;
    adj_list.v[v].etotal++;
    adj_list.etotal++;
}

/*
 * Linked Forward List
 * @type Data Structure
 */
struct Adjedge {
    int v, pre;
    WTYPE w;
};
Adjedge adj_edge[2 * MAXE];
int head[MAXV];
int h_head(0);
// add u->v
void add_edge(int u, int v, WTYPE w) {
    adj_edge[++h_head].v = v;
    adj_edge[h_head].w = w;
    adj_edge[h_head].pre = head[u];
    head[u] = h_head;
}
#define tranverse_edge_of_u(i, x) \
    for (int i = head[(x)]; i; i = adj_edge[i].pre)

/**
 * Disjoint Set
 * @type Data Structure
 */
int rnk[MAXV];  // Merge with rank
int DJS[MAXV];  // Disjoint set
int num[MAXV];  // Attr: num, max, min, rel(or use factDJS)

static int find(int x) {
    if (DJS[x] == x)
        return x;
    else
        return (DJS[x] = find(DJS[x]));
}

int eq(int x, int y) { return (find(x) == find(y)); }

// if weighted-DJS is used, do not use rank optimization
void merge(int x, int y) {
    x = find(x);
    y = find(y);
    if (x == y) return;
    if (rnk[x] < rnk[y]) {
        DJS[x] = y;
        num[y] = num[x] + num[y];  // Attr alters
    } else {
                DJS[y] = x;
        num[x] = num[y] + num[x];  // Attr alters
    }
    if (rnk[x] == rnk[y]) rnk[x]++;
}

void DJS_initial(int n) {
    for (int i = 1; i <= n; ++i) {
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
 * @return weight of minimal spinning tree (-1 if not exist)
 * @timecomplex O(Elog(E))
 */
WTYPE Kruskal(int n, int m) {
    WTYPE ret = 0;
    int cnt = 0;
    DJS_initial(n);
    sort(edge_set + 1, edge_set + m + 1);
    for (int i = 1; i <= m; ++i)
        if (!eq(edge_set[i].u, edge_set[i].v)) {
            cnt++;
            merge(edge_set[i].u, edge_set[i].v);
            ret += edge_set[i].w;
        }

    if (cnt == n - 1)
        return ret;
    else
        return -1;
}

/**
 * Prim
 * @type Algorithm for minimal spinning tree
 * @param 1.adjacency list 2.adj_mtrx 3. Fibonacci Heap
 * @param n vertices number
 * @return weight of minimal spinning tree (-1 if not exist)
 * @timecomplex 1.O(Elog(V)) 2.O(V^2) 3. O(E+VlgV)
 */
WTYPE mincost[MAXV];
priority_queue<pTi, vector<pTi>, greater<pTi> > q;
bool vis[MAXV];
WTYPE Prim(int n) {
    WTYPE ret = 0;
    int cnt = 0;
    memset(mincost,0x3f,sizeof(WTYPE)*(n+1));
    memset(vis,0,sizeof(bool)*(n+1));
    while (!q.empty()) q.pop();     
    mincost[1] = 0;
    q.push(make_pair(0, 1));

    while (!q.empty() && cnt < n) {
        while (vis[q.top().second]) q.pop();
        WTYPE w = q.top().first;
        int u = q.top().second;
        q.pop();
        cnt++;
        ret += w;
        vis[u] = true;
        tranverse_edge_of_u(i, u) {
            if (adj_edge[i].w < mincost[adj_edge[i].v]) {
                mincost[adj_edge[i].v] = adj_edge[i].w;
                q.push(make_pair(adj_edge[i].w, adj_edge[i].v));
            }
        }
    }

    if (cnt == n)
        return ret;
    else
        return -1;
}


/**
 * Miniaml Spanning Tree Related
 * @type Algorithm List
 * -Identify the uniqueness of MST
 */

/*Identify the uniqueness of MST
1.Kruskal and mark the chosen edges
2.Sort with second keyword that marked one are sorted behind
3.Kruskal and mark the chosen edges with another mark
4.Jugde whether marks of all edges are same respectively
*/


/**
 * Cayley Tree Formula
 * @type Formula to count tree of order n
 * @param n order of graph
 * @return number of different tree
 */
int Cayley(int n) { return fp(n, n - 2); }


/**
 * Matrix Tree Theorem
 * @type Algorithm to count tree of any graph
 * @param n
 * @param deg
 * @param adj_mtrx
 * @return number of different tree
 */
double C_mtrx[MAXV][MAXV];
int number_of_spinning_tree(int n) {
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            if (i != j) {
                C_mtrx[i][j] = -adj_mtrx[i][j];
            } else {
                C_mtrx[i][j] = deg[i];
            }
        }
    }
    return (int)(det(n - 1, C_mtrx) + 0.5);
}


/**
 * DFS
 * @type Algorithm
 * @param Linked Forward List
 * @return DFS_pre father of node in DFS tree
 * @return ordL,ordR DFS order of node in DFS tree
 * @return b,f timestack of beginning and finish time
 * @return cc connection component
 * @return simple_way ways to a node from root in DFS tree
 * XXX: UNCHECKED
 */
int dfs_time,ord_time,counter;
bool DFS_vis[MAXV];
int DFS_pre[MAXV],ordL[MAXV],ordR[MAXV],b[MAXV],f[MAXV];
bool connected;
//Only for undirected-graph
int cc[MAXV];
bool is_dag;
//Only for DAG
int simple_way[MAXV];

void DFS_visit(int u) {
    cc[u] = counter;
    DFS_vis[u] = true;
    ord_time++; dfs_time++;
    b[u] = dfs_time;
    ordL[u] = ord_time;

    tranverse_edge_of_u(i,u){
        if (vis[adj_edge[i].v]){
            if (adj_edge[i].v!=u) is_dag = false;
            continue; //if a tree is given, pre==u can be used
        }
        simple_way[adj_edge[i].v] += simple_way[u];
        DFS_vis[adj_edge[i].v] = u;
        DFS_visit(adj_edge[i].v);
    }

    dfs_time++;
    f[u] = dfs_time;
    ordR[u] = ord_time;
    
}

void DFS(int n){
    is_dag = true;
    for (int i=0;i<n;++i){
        DFS_vis[i] = false;
        DFS_pre[i] = 0;
    }
    int cnt = 0;
    ord_time = dfs_time = counter = 0;
    for (int i=0;i<n;++i){
        if (!DFS_vis[i]){
            cnt++;
            cc[i] = ++counter;
            simple_way[i] = 1;
            DFS_visit(i);
        }
    }
    connected = cnt>1?false:true;
}



/**
 * Topological sort
 * @type Algorithm
 * @param graph
 * @return topo_sort topological order of vertices
 */
int topo_sort[MAXV];
bool topocmp(int a, int b){
    return f[a]>f[b];
}
void Topological_sort(int n){
    DFS(n);
    sort(topo_sort, topo_sort+n, topocmp);
}


/**
 * Strongly Connected Component
 * @WTYPE Algorithm Tarjan, Kosaraju is deprecated
 * @param graph
 * @return scc strong connected component
 * @return new_dag Linked Forward List for new graph
 */

priority_queue<pTi, vector<pTi>, greater<pTi> > Q;
WTYPE dis[MAXV];
//bool vis[MAXV];
int Dijkstra(int S, int T){
    while (!Q.empty()) Q.pop();
    memset(dis, 0x3f, sizeof(dis));
    memset(vis, 0, sizeof(vis));
    Q.push(make_pair(0, S));
    dis[S] = 0;

    while (!Q.empty()){
        pTi cur = Q.top();
        Q.pop();
        if (dis[cur.second] < cur.first) continue;
        for (int i = head[cur.second]; i; i = adj_edge[i].pre)
            if (dis[cur.second] + adj_edge[i].w < dis[adj_edge[i].v])
                Q.push(make_pair(dis[adj_edge[i].v]=dis[cur.second]+adj_edge[i].w, adj_edge[i].v));
    }

    return dis[T];
}


#endif

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
typedef pair<int, int> Pii;
int n,m;


/**
 * Edge Set
 */
struct Edge {
    int u, v;
    WTYPE w;
    bool operator<(const Edge &t) const { return w < t.w; }
};
Edge edge_set[MAXE];  // start index from 1
int edge_h;


/**
 * Adjacency Matrix
 */
WTYPE adjmtrx[MAXV][MAXV];
int deg[MAXV];


/**
 * Adjacency List
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
void addedge_adjlist(int u, int v, WTYPE w) {
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
 */
struct Adjedge {
    int v, pre;
    WTYPE w;
};
Adjedge adjedge[2 * MAXE];
int head[MAXV];
int head_h(0);
// add u->v
void addedge(int u, int v, WTYPE w) {
    adjedge[++head_h].v = v;
    adjedge[head_h].w = w;
    adjedge[head_h].pre = head[u];
    head[u] = head_h;
}
void init_adjedge(){
    head_h = 0;
    memset(adjedge,-1,sizeof(adjedge));
}
    

/**
 * Disjoint Set
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
        for (int i = head[u]; i; i = adjedge[i].pre) {
            if (adjedge[i].w < mincost[adjedge[i].v]) {
                mincost[adjedge[i].v] = adjedge[i].w;
                q.push(make_pair(adjedge[i].w, adjedge[i].v));
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
                C_mtrx[i][j] = -adjmtrx[i][j];
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

    for (int i = head[u]; i; i = adjedge[i].pre){
        if (vis[adjedge[i].v]){
            if (adjedge[i].v!=u) is_dag = false;
            continue; //if a tree is given, pre==u can be used
        }
        simple_way[adjedge[i].v] += simple_way[u];
        DFS_vis[adjedge[i].v] = u;
        DFS_visit(adjedge[i].v);
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
 * Dijkstra
 * @type Algorithm for SSSP
 * param S: source
 * param T: destination
 */

priority_queue<pTi, vector<pTi>, greater<pTi> > Q;
WTYPE dis[MAXV];
int Dijkstra(int S, int T){
    while (!Q.empty()) Q.pop();
    memset(dis, 0x3f, sizeof(dis));
    Q.push(make_pair(0, S));
    dis[S] = 0;

    while (!Q.empty()){
        pTi cur = Q.top();
        Q.pop();
        if (dis[cur.second] < cur.first) continue;
        for (int i = head[cur.second]; i; i = adjedge[i].pre)
            if (dis[cur.second] + adjedge[i].w < dis[adjedge[i].v])
                Q.push(make_pair(dis[adjedge[i].v]=dis[cur.second]+adjedge[i].w, adjedge[i].v));
    }

    return dis[T];
}


/**Tarjan
 * @type Algorithm on undigraph for BCC, cut-vertex and bridge
 */

int low[MAXV],dfn[MAXV],tindex;
stack<int> S;

int iscut[MAXV];
pTi bridge[MAXE];
int bridge_h,number_of_cut,number_of_bridge,number_of_bcc;
int bccno[MAXV];
// For bcc, the cut-vertice may not be in every bcc
void tarjan_dfs(int cur, int pre) {
    dfn[cur] = low[cur] = ++tindex;
    int child = 0;
    S.push(cur);
    for (int i = head[cur]; i; i = adjedge[i].pre){
        int v = adjedge[i].v;
        if (!dfn[v]) {
            tarjan_dfs(v, cur);
            low[cur] = min(low[cur],low[v]);

            child ++;
            if (low[v]>=dfn[cur]){
                if (!iscut[cur]) number_of_cut ++;
                iscut[cur] ++;
            }

            if (low[v]>dfn[cur]){
                number_of_bridge ++;
                bridge[bridge_h].first = min(v,cur);
                bridge[bridge_h++].second = max(v,cur);
            }
        } else if (dfn[v] < dfn[cur] && v!=pre) {
            low[cur] = min(low[cur], dfn[v]);
        }
    }

    if (low[cur] == dfn[cur]){
        number_of_bcc++;
        int v;
        do{
            bccno[v = S.top()] = number_of_bcc;
            S.pop();
        } while (v!=cur);
    }

    if (pre == -1) {
        if (child == 1) {
            number_of_cut--;
            iscut[cur] = false;
        } else {
            iscut[cur] --;
        }
    }
}

void tarjan_undigraph(){
    tindex = bridge_h = number_of_cut = number_of_bridge = 0;
    for (int i=0;i<n;++i){
        if (!dfn[i]){
            tarjan_dfs(i, -1);
        }
    }
}


/** Tarjan SCC
 * @type Algorithm on digraph for SCC
 */
int sccno[MAXV],scccnt[MAXV];
int number_of_SCC;
bool instack[MAXV];
void SCC_dfs(int cur, int pre) {
    dfn[cur] = low[cur] = ++tindex;
    instack[cur] = true;
    S.push(cur);
    for (int i = head[cur]; i; i = adjedge[i].pre){
        int v = adjedge[i].v;
        if (!dfn[v]) {
            tarjan_dfs(v, cur);
            low[cur] = min(low[cur],low[v]);
        } else if (instack[v]) {
            low[cur] = min(low[cur],dfn[v]);
        }
    }

    if (dfn[cur] == low[cur]){
        number_of_SCC++;
        int j;
        do {
            num[number_of_SCC] ++;
            sccno[j = S.top()] = number_of_SCC;
            S.pop();
            instack[j] = false;
        } while (j!=cur);
    }
}

void tarjan_SCC(){
    number_of_SCC = 0;
    for (int i=1;i<=n;++i){
        if (!dfn[i]){
            SCC_dfs(i, -1);
        }
    }
}


#endif

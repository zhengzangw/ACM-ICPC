#include <iostream>
#include <cstdio>
#include <vector>
#include <cmath>
const int maxn=100;
const int INF =1<<20;
int n,m,t,ans;
struct edge{int to,length;};
using namespace std;
vector<edge> G[maxn];
bool centroid[maxn];
int subtree_size[maxn];
int k;

int compute_subtree_size(int v,int p)
{
    int c = 1;
    for (int i = 0;i<G[v].size();i++){
        int w = G[v][i].to;
        if (w == p || centroid[w]) continue;
        c += compute_subtree_size(w,v);
    }
    subtree_size[v] = c;
    return c;
}

pair<int,int> search_controid(int v,int p,int t)
{
    pair<int,int> res = make_pair(INF,-1);
    int s=1, m=0;
    for (int i=0;i<G[v].size();i++){
        int w = G[v][i].to;
        if (w==p || centroid[w]) continue;
        res = min(res,search_controid(w,v,t));
        m = max(m,subtree_size[w]);
        s+=subtree_size[w];
    }
    m = max(m,t-s);
    res = min(res,make_pair(m,v));
    return res;
}

void enumerate_paths(int v, int p, int d, vector<int> &ds)
{
    ds.push_back(d);
    for (int i=0;i<G[v].size();i++){
        int w=G[v][i].to;
        if (w==p || centroid[w]) continue;
        enumerate_paths(w,v,d+G[v][i].length, ds);
    }
}

int count_pair(vector<int> &ds)
{
    int res=0;
    sort(ds.begin(),ds.end());
    int j=(int)ds.size();
    for (int i=0;i<ds.size();i++){
        while (j>0 && ds[i]+ds[j-1]>k ) --j;
        res += j-(j>i?1:0);
    }
    return res/2;
}

void solve_subproblem(int v)
{
    compute_subtree_size(v,-1);
    int s = search_controid(v,-1,subtree_size[v]).second;
    centroid[s]=true;
    
    for (int i=0;i<G[s].size();i++){
        if (centroid[G[s][i].to]) continue;
        solve_subproblem(G[s][i].to);
    }
    
    vector<int> ds;
    ds.push_back(0);
    for (int i=0;i<G[s].size();i++){
        if (centroid[G[s][i].to])continue;
        vector<int> tds;
        enumerate_paths(G[s][i].to,s,G[s][i].length,tds);
        ans-=count_pair(tds);
        ds.insert(ds.end(),tds.begin(),tds.end());
    }
    ans += count_pair(ds);
    centroid[s]=false;
}

void solve()
{
    ans = 0;
    solve_subproblem(1);
    printf("%d\n",ans );
}

int main(){
    int n;
    cin >> n >> k;
    int a,b,l;
    for (int i=1;i<n;i++){
        cin >> a >> b >> l;
        edge tm; tm.to = b; tm.length = l;
        G[a].push_back(tm);
        tm.to = a; tm.length = l;
        G[b].push_back(tm);
    }
    memset(centroid,0,sizeof(centroid));
    solve();
    return 0;
}

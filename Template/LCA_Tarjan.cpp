#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N = 1 << 17;
const int INF = 1 << 27;
int n,m,t,x,y,z;
struct Edge
{
    int u,v,c,next;
};
Edge edge[MAX_N];
int p[MAX_N],h;
bool vis[MAX_N];
int st[MAX_N][32],f[MAX_N][32],deep[MAX_N];
struct que{
    int ind,u;
};
que q[MAX_N][MAX_N];
int le[MAX_N];
int par[MAX_N],ans[MAX_N],ran[MAX_N];

void addedge(int u,int v,int c){
    edge[h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h++;
}

int find(int x){
    if (par[x] == x){
        return x;
    } else {
        return (par[x] = find(par[x]));
    }
}

void Tarjan(int u){
    
    for (int i=1;i<=le[u];i++)
        if (vis[q[u][i].u])
            ans[q[u][i].ind]=find(q[u][i].u);
    
    vis[u]=true;
    
    for (int i=p[u];i!=-1;i=edge[i].next)
        if (!vis[edge[i].v]) {
            Tarjan(edge[i].v);
            par[edge[i].v]=u;
            find(edge[i].v);
        }
    
}
/*If we need to calculate the length of the path, then
改造并查集,定义dis[i]数组,保存i到getf(i)的距离
定义Deep[i]数组,表示i节点的深度,DFS时顺便更新depp[i];
定义Sum[I]数组,表示从根节点到I深度节点的距离.因为在LCA Tarjan算法中 ,LCA(设为X) 必然在DFS路径上,所以X到I的距离为sum[deep[I]]-sum[Deep[X]]
响应时,返回值为：dis[A]-sum[deep[find(A)]]+sum[Deep[B]];*/
int main(){
    cin >> n;
    h=0;
    fill(p,p+n+1,-1);
    for (int i=0;i<n-1;i++){
        cin >> x >> y >> z;
        addedge(x,y,z);
        addedge(y,x,z);
    }
    
    for (int i=1;i<=n;i++) par[i]=i;
    
    cin >> m;
    for (int i=1;i<=m;i++){
        cin >> x >> y;
        q[x][++le[x]].u=y;
        q[x][le[x]].ind=i;
        q[y][++le[y]].u=x;
        q[y][le[y]].ind=i;
        
    }
    cout << "DONE!\n";
    fill(vis,vis+n+1,false);
    Tarjan(1);
    for (int i=1;i<=m;i++) cout<<ans[i]<<'\n';
    return 0;
}
#include <iostream>
#include <cstdio>
#include <cmath>
using namespace std;
const int MAX_N = 1 << 10;
const int INF = 1 << 27;
int n,m,h,s,t;
int p[MAX_N];
bool used[MAX_N];
struct Edge{
    int u,v,c,next;
};
Edge edge[MAX_N*3];

void addedge(int u, int v,int c)
{
    edge[h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h++;
    edge[h].u = v; edge[h].v = u; edge[h].c = 0;
    edge[h].next = p[v]; p[v] = h++;
}

void init()
{
    int T;
    scanf("%d",&T);
    while (T--){
        h = 0;
        scanf("%d%d",&n,&m);
        scanf("%d%d",&s,&t);
        fill(p,p+n+1,-1);
        for (int i=0;i<m;i++){
            int a,b,c;
            scanf("%d%d%d",&a,&b,&c);
            addedge(a,b,c);
        }
    }
}

int dfs(int v, int t, int f)
{
    if (v==t) return f;
    used[v] = true;
    for (int i = p[v];i!=-1;i=edge[i].next){
        int e = edge[i].v, w = edge[i].c;
        if (!used[e] && w>0)
        {
            int d=dfs(e,t,min(f,w));
            if (d>0) {
                edge[i].c -=d;
                edge[i^1].c +=d;
                return d;
            }
        }
    }
    return 0;
}

int Max_Flow(int s, int t)
{
    int flow=0;
    while(true){
        memset(used,0,sizeof(used));
        int f = dfs(s,t,INF);
        if (f == 0) return flow;
        flow += f;
    }
}

int main(){
    init();
    printf("%d\n",Max_Flow(s,t));
    return 0;
}
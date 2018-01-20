#include<cstdio>
#include<cstring>
#include<queue>
#include<cmath>

using namespace std;
const int Maxn = 1100;
const int Maxm = 12000;
const int INF = 1 << 26;
struct Edge{
    int u,v,c,next;
};
Edge edge[Maxm];
int n,m,i,h,s,t,T;
int dep[Maxn],p[Maxn];
queue <int> q;

void addedge(int u,int v,int c){
    edge[h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h++;
    
    edge[h].u = v; edge[h].v = u; edge[h].c = 0;
    edge[h].next = p[v]; p[v] = h++;
}

int bfs(int s,int t){
    while (!q.empty()) q.pop();
    memset(dep,-1,sizeof(dep));
    q.push(s); dep[s] = 0;
    while (!q.empty()){
        int u = q.front(); q.pop();
        for (int i = p[u];i!=-1;i=edge[i].next){
            int v = edge[i].v;
            if (dep[v] == -1 && edge[i].c>0)
            {
                dep[v] = dep[u] + 1;
                q.push(v);
            }
        }
    }
    return dep[t] != -1;
}

int dfs(int s,int t,int m){
    int r=0;
    if (s == t) return m;
    for (int i = p[s];i!=-1 && r < m;i = edge[i].next){
        int v = edge[i].v;
        if (edge[i].c > 0 && dep[v] == dep[s]+1){
            int x = min(edge[i].c,m-r);
            x = dfs(v,t,x);
            r += x;
            edge[i].c -= x;
            edge[i^1].c += x;
        }
    }
    if (!r) dep[s] = -2;
    return r;
}

int dinic(int s,int t){//O(n2m)
    int maxflow = 0,ex_f;
    while (bfs(s,t)){
        while (ex_f = dfs(s,t,INF))
            maxflow += ex_f;
    }
    return maxflow;
}
int main(){
    scanf("%d",&T);
    while (T--){
        memset(p,-1,sizeof(p));
        h = 0;
        scanf("%d%d",&n,&m);
        scanf("%d%d",&s,&t);
        for (int i=0;i<m;i++){
            int a,b,c;
            scanf("%d%d%d",&a,&b,&c);
            //Calculate the shortest path 
            //addedge(a,b,c*337+1);
            addedge(a,b,c);
        }
        
        //printf("%d\n",dinic(s,t)%337);
        printf("%d\n",dinic(s,t));
    }
    return 0;
}

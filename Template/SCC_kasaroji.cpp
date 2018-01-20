#include<cstdio>
#include<cstring>
#include<queue>
#include<cmath>

using namespace std;
const int Maxn = 1100;
const int Maxm = 12000;
struct Edge{
    int u,v,c,next;
};
Edge edge[Maxm],redge[Maxm];
int n,m,i,h,s,t,T,k,v,u,rh;
int p[Maxn];
bool instack[Maxn];
int used[Maxn],se[Maxn],low[Maxn],rp[Maxn];
int top,ind;

void addedge(int u,int v,int c){
    edge[++h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h;
}

void r_addedge(int u,int v,int c){
    redge[++rh].u = u; redge[rh].v = v; redge[rh].c = c;
    redge[rh].next = rp[u]; rp[u] = rh;
}

void dfs(int v){
    used[v]=true;
    for (int i=p[v];i!=-1;i=edge[i].next)
        if (!used[edge[i].v]) dfs(edge[i].v);
    se[++top]=v;
}

void rdfs(int v,int k){
    used[v]=true;
    low[v]=k;
    for (int i=rp[v];i!=-1;i=redge[i].next)
        if (!used[redge[i].v]) rdfs(redge[i].v,k);
}

int solve()
{
    memset(used,0,sizeof(used));
    top=0;
    for (int v=1;v<=n;v++)
        if (!used[v]) dfs(v);
    memset(used,0,sizeof(used));
    int k=0;
    for (int v=top;v>0;v--)
        if (!used[se[v]]) rdfs(se[v],k++);
    return k;
}

int main(){
    h=0; rh=0;
    fill(p,p+Maxn,-1);
    fill(rp,rp+Maxn,-1);
    scanf("%d%d",&n,&m);
    for (int i=0;i<m;i++){
        scanf("%d%d",&v,&u);
        addedge(v,u,1);
        r_addedge(u,v,1);
    }
    
    printf("%d\n",solve());
    return 0;
}

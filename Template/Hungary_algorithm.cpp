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
int n,m,i,h,s,t,T,k,v,u;
int dep[Maxn],p[Maxn],from[Maxn];
bool use[Maxn];

void addedge(int u,int v,int c){
    edge[++h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h;
}

bool match(int x)
{
    for (int i=p[x];i!=-1;i=edge[i].next){
        if (!use[edge[i].v]){
            use[edge[i].v]=true;
            if (from[edge[i].v] == -1||match(from[edge[i].v])){
                from[edge[i].v] = x;
                return true;
            }
        }
    }
    return false;
}

int Hun(){//O(VE)
    int tot=0;
    fill(from,from+m+1,-1);
    for (int i=1;i<=m;++i){
        memset(use,0,sizeof(use));
        if (match(i)) tot++;
    }
    return tot;
}

int main(){
    scanf("%d",&T);
    h=0;
    fill(p,p+Maxn,-1);
    while (T--){
        scanf("%d%d",&m,&k);
        for (int i=0;i<k;i++){
            scanf("%d%d",&v,&u);
            addedge(v,u,1);
        }
    }
    printf("%d\n",Hun());
    return 0;
}

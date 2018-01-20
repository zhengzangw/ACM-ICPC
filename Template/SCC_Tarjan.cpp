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
bool instack[Maxn];
int DFN[Maxn],LOW[Maxn],pre[Maxn],Stap[Maxn],Belong[Maxn];
int Bcnt,Stop,ind;

void addedge(int u,int v,int c){
    edge[++h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h;
}

void tarjan(int i)
{
    int j;
    DFN[i]=LOW[i]=++ind;
    instack[i]=true;
    Stap[++Stop]=i;
    for (int v=p[i];v!=-1;v=edge[v].next){
        j=edge[v].v;
        if(!DFN[j]){
            tarjan(j);
            if(LOW[j]<LOW[i]) LOW[i]=LOW[j];
        } else if (instack[j] && DFN[j]<LOW[i])
            LOW[i]=DFN[j];
    }
    
    if(DFN[i]==LOW[i]){
        Bcnt++;
        do{
            j=Stap[Stop--];
            instack[i]=false;
            Belong[j]=Bcnt;
        }
        while (j!=i);
    }
}

void solve()
{
    int i;
    Stop=Bcnt=ind=0;
    memset(DFN,0,sizeof(DFN));
    for (i=1;i<=n;i++)
        if (!DFN[i]) tarjan(i);
}

int main(){
    h=0;
    fill(p,p+Maxn,-1);
    scanf("%d%d",&n,&m);
    for (int i=0;i<m;i++){
        scanf("%d%d",&v,&u);
        addedge(v,u,1);
    }
    
    solve();
    printf("%d\n",Bcnt);
    return 0;
}

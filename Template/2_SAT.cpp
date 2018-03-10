#include<cstdio>
#include<cstring>
#include<cmath>
#include<iostream>

using namespace std;
const int Maxn = 1100;
const int Maxm = 12000;
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
    fill(DFN,DFN+2*n+1,0);
    for (i=1;i<=2*n;i++)
        if (!DFN[i]) tarjan(i);
}

void judge()
{
    for (int i=1;i<=n;i++)
        if (Belong[i]==Belong[i+n]){
            printf("NO");
            return;
        }
    
    printf("YES\n");
    /*for (int i=1;i<=n;i++)
        if (Belong[i]>Belong[i+n]) printf("%d:true\n",i);
        else printf("%d:false\n",i);*/
    //There are some problems
}

int chu(int x)
{
    if (x<=n) return x+n;
    else return x-n;
}

int main(){
    h=0;
    scanf("%d%d",&n,&m);
    fill(p,p+n*2+1,-1);
    for (int i=0;i<m;i++){
        scanf("%d%d",&v,&u);
        if (v<0) v=-v+n;
        if (u<0) u=-u+n;
        addedge(chu(v),u,0);
        addedge(chu(u),v,0);
    }
    
    solve();
    judge();
    
    return 0;
}

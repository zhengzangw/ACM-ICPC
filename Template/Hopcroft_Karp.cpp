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
int dep[Maxn],p[Maxn],from[Maxn],mx[Maxn],my[Maxn],dx[Maxn],dy[Maxn];
bool use[Maxn];
queue<int> que;

void addedge(int u,int v,int c){
    edge[++h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h;
}

bool find(int x)
{
    for (int i=p[x];i!=-1;i=edge[i].next){
        if (!use[edge[i].v] && dy[edge[i].v] == dx[x] + 1){
            use[edge[i].v]=true;
            if (!my[edge[i].v]||find(my[edge[i].v])){
                my[edge[i].v] = x;
                mx[x] = edge[i].v;
                return true;
            }
        }
    }
    return false;
}

int matching(){//O(VE^0.5)
    memset(mx,0,sizeof(mx));
    memset(my,0,sizeof(my));
    int ans=0;
    while (true){
        bool flag=false;
        while (!que.empty()) que.pop();
        memset(dx,0,sizeof(dx));
        memset(dy,0,sizeof(dy));
        for (int i=1; i<=n; i++)
            if (!mx[i]) que.push(i);
        while (!que.empty()) {
            int u=que.front(); que.pop();
            for (int i=p[u];i!=-1;i=edge[i].next)
                if (!dy[edge[i].v]){
                    dy[edge[i].v] = dx[u] + 1;
                    if (my[edge[i].v]){
                        dx[my[edge[i].v]] = dy[edge[i].v] + 1;
                        que.push(my[edge[i].v]);
                    } else flag = true;
                }
        }
        
        if (!flag) break;
        memset(use,0,sizeof(use));
        for(int i=1; i<=n;i++)
            if (!mx[i] && find(i)) ans++;
                }
    return ans;
}

int main(){
    scanf("%d",&T);
    h=0;
    fill(p,p+Maxn,-1);
    while (T--){
        scanf("%d%d",&n,&k);
        for (int i=0;i<k;i++){
            scanf("%d%d",&v,&u);
            addedge(v,u,1);
        }
    }
    printf("DONE!\n");
    printf("%d\n",matching());
    return 0;
}

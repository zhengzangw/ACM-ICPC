#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N = 1 << 17;
int n,m,t,x,y,z;
struct Edge
{
    int u,v,c,next;
};
Edge edge[MAX_N];
int p[MAX_N],h,ans;
bool vis[MAX_N];
int st[MAX_N][32],f[MAX_N][32],g[MAX_N][32],deep[MAX_N];

void addedge(int u,int v,int c){
    edge[h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h++;
}

void dfs(int a,int dep)
{
    vis[a] = true;
    int i=1,t;
    deep[a] = dep;
    while ((1<<i) <= dep){
        f[a][i]=f[f[a][i-1]][i-1];
        g[a][i]=f[f[a][i-1]][i-1]+g[a][i-1];
        i++;
    }
    for (i=p[a];i!=-1;i=edge[i].next){
        t = edge[i].v;
        if (!vis[t]) {
            f[t][0]=a;
            g[t][0]=edge[i].c;
            dfs(t,dep+1);
        }
    }
}

int jump(int &a,int dep)
{
    int i=0,ret=0;
    while (dep>0){
        if (dep & 1) {
            ret += g[a][i];
            a= f[a][i];
        }
        i++;
        dep >>=1;
    }
    return ret;
}

int query(int a,int b)
{
    int ret=0;
    if (deep[a]>deep[b]) ret+= jump(a,deep[a]-deep[b]);
    else ret+=jump(b,deep[b]-deep[a]);
    
    if (a==b) {ans=a; return ret;}
    
    for(int i=log2(n);i>=0;i--){
        if (f[a][i]==0) continue;
        if (f[a][i]!=f[b][i]){
            ret += g[a][i];
            ret += g[b][i];
            a = f[a][i];
            b = f[b][i];
        }
    }
    
    if (a!=b) ret += g[a][0]+g[b][0], a=f[a][0];
    
    ans = a;
    return ret;
}

int main(){
	//initial
    cin >> n;
    h=0;
    fill(p,p+n+1,-1);
    for (int i=0;i<n-1;i++){
        cin >> x >> y >> z;
        addedge(x,y,z);
        addedge(y,x,z);
    }
    
    dfs(1,0);
    
    //test
    cin >> m;
    while(m--){
        cin >> x >> y;
        cout << query(x,y) << " " << ans << "\n";
    }

    return 0;
}
#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N = 100;
int n,m,t,x,y,z;
struct Edge
{
    int u,v,c,next;
};
Edge edge[MAX_N];
int p[MAX_N],h,ans;
bool vis[MAX_N];
int st[MAX_N][32],f[MAX_N][32],g[MAX_N][32],deep[MAX_N];
int id[MAX_N],vs[MAX_N],depth[MAX_N];
int smin[MAX_N][MAX_N],sind[MAX_N][MAX_N];

void addedge(int u,int v,int c){
    edge[h].u = u; edge[h].v = v; edge[h].c = c;
    edge[h].next = p[u]; p[u] = h++;
}

void dfs(int v,int q, int d, int &k)
{
    id[v] = k; //first time visit v
    vs[k] = v; //visiting sequence
    depth[k++] = d;
    for (int i=p[v]; i!=-1; i=edge[i].next){
        if (edge[i].v != q){
            dfs(edge[i].v,v,d+1,k);
            vs[k]=v;
            depth[k++] = d;
        }
    }
}

void rmq_init(int n){
    for (int i=1;i<=n;i++) {smin[i][0]=depth[i]; sind[i][0]=i;}
    int m = ceil(log2(n));
    for (int j=1; j<m;j++)
        for (int i=0; i+ (1<<j)-1 < n;i++)
            if (smin[i][j-1]>smin[i+(1<<(j-1))][j-1]){
                smin[i][j]=smin[i+(1<<(j-1))][j-1];
                sind[i][j]=sind[i+(1<<(j-1))][j-1];
            } else {
                smin[i][j]=smin[i][j-1];
                sind[i][j]=sind[i][j-1];
            }
}

int query(int L,int R){
    int k = log2(R-L+1);
    if (smin[L][k]>smin[R-(1 << k) + 1][k]) return sind[R-(1<<k)+1][k];
    else return sind[L][k];
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
    
    int k=1;
    dfs(1,-1,0,k);
    rmq_init(n*2);
    //test
    cin >> m;
    while(m--){
        cin >> x >> y;
        int xx=min(id[x],id[y]);
        int yy=max(id[x],id[y]);
        cout << vs[query(xx,yy)] << "\n";
    }
    
    return 0;
}
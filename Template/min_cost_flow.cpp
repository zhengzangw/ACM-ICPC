#include <iostream>
#include <cstring>
#include <cmath>
#include <queue>
using namespace std;
typedef pair<int,int> P;
const int MAX_N = 1 << 17;
const int INF = 1 << 27;
int n,m,t,h;

int dist[MAX_N],prevv[MAX_N],preve[MAX_N],hs[MAX_N];

int p[MAX_N];
struct Edge{
    int u,v,cap,cost,next;
};
Edge edge[MAX_N*2];
void addedge(int u,int v,int cap,int cost)
{
    edge[h].u=u; edge[h].v=v; edge[h].cap=cap;
    edge[h].cost=cost; edge[h].next=p[u]; p[u]=h++;
    
    edge[h].u=v; edge[h].v=u; edge[h].cap=0;
    edge[h].cost=-cost; edge[h].next=p[v]; p[v]=h++;
}

/*int min_cost_flow_B(int s, int t, int f){ //Based On Bellman-Ford O(FVE)
 int res = 0;
 while (f>0){
 fill(dist,dist+n+1,INF);
 dist[s]=0;
 bool update = true;
 while (update) {
 update = false;
 for (int v=1; v<=n; v++){
 if (dist[v] == INF) continue;
 for (int i=p[v];i!=-1;i=edge[i].next){
 Edge e=edge[i];
 if (e.cap>0 && dist[e.v]>dist[v]+e.cost){
 dist[e.v] = dist[v] + e.cost;
 prevv[e.v]=v;
 preve[e.v]=i;
 update = true;
 }
 }
 }
 }
 
 
 if (dist[t] == INF) return -1;
 
 int d=f;
 for (int v = t; v!=s; v=prevv[v]){
 d = min(d,edge[preve[v]].cap);
 }
 f -= d;
 res += d*dist[t];
 for (int v = t; v!=s; v=prevv[v]){
 edge[preve[v]].cap -= d;
 edge[preve[v]^1].cap += d;
 }
 }
 return res;
 }*/

int min_cost_flow(int s, int t, int f){ //Based On Dijskra O(FV2)
    int res = 0;
    fill(hs,hs+n+1,0);
    
    while (f>0){
        priority_queue<P,vector<P>,greater<P> > que;
        fill(dist,dist+n+1,INF);
        dist[s]=0;
        que.push(P(0,s));
        while (!que.empty()){
            P pp = que.top(); que.pop();
            int v = pp.second;
            if (dist[v] < pp.first) continue;
            for (int i=p[v];i!=-1;i=edge[i].next){
                Edge e=edge[i];
                if (e.cap>0 && dist[e.v] > dist[v]+e.cost+hs[v]-hs[e.v]){
                    dist[e.v]=dist[v]+e.cost+hs[v]-hs[e.v];
                    prevv[e.v]=v;
                    preve[e.v]=i;
                    que.push(P(dist[e.v],e.v));
                }
            }
        }
        
        
        if (dist[t] == INF) return -1;
        for (int v=1;v<=n;v++) hs[v] += dist[v];
        
        int d=f;
        for (int v = t; v!=s; v=prevv[v]){
            d = min(d,edge[preve[v]].cap);
        }
        f -= d;
        res += d*hs[t];
        for (int v = t; v!=s; v=prevv[v]){
            edge[preve[v]].cap -= d;
            edge[preve[v]^1].cap += d;
        }
    }
    return res;
}

int main(){
    
    int cap,cost,u,v,s,t,f;
    cin >> n >> m;
    cin >> s >> t >> f;
    h = 0;
    fill(p,p+n+1,-1);
    for (int i=0;i<m;i++){
        cin >> u >> v >> cap >> cost;
        addedge(u,v,cap,cost);
    }
    
    cout << "The Answer is: " << min_cost_flow(s,t,f) << '\n';
    return 0;
}
#include <iostream>
#include <string>
#include <vector>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,f,V,E,t;
int cost[MAX_N][MAX_N];
int mincost[MAX_N],used[MAX_N];

int MST(){ //O(E*log(V))
    for (int i=1;i<=V;++i){
        mincost[i] =INF;
        used[i] = false;
    }
    mincost[1] = 0;
    int res = 0;
    
    while (true){
        int v=-1;
        for (int u=1;u<=V;u++)
            if (!used[u] && (v==-1 || mincost[u]<mincost[v])) v = u;
        
        if (v==-1) break;
        used[v] = true;
        res += mincost[v];
        
        for (int u=1; u<=V; u++){
            mincost[u] = min(mincost[u],cost[v][u]);
        }
    }
    
    return res;
}

void init(){
    cin >> V >> E;
    for (int i=1;i<=V;i++){
        for (int j=1;j<=V;j++){
            cost[i][j]=INF;
            cost[j][i]=INF;
        }
    }
    for (int i=1;i<=V;i++) cost[i][i]=0;
    
    for (int i=0;i<E;i++){
        cin >> v >> u >> w;
        cost[v][u] = w;
        cost[u][v] = w;
    }
    
}

int main(){
    init();
    cout<<"DONE!\n";
    cout << MST();
    return 0;
}
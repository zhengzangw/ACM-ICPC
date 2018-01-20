#include <iostream>
#include <string>
#include <vector>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,f,V,E,t;
int cost[MAX_N][MAX_N];

void shortest_path(){ //O(V^3)
    for (int k=1;k<=V;k++)
        for (int i=1;i<=V;i++)
            for (int j=1;j<=V;j++)
                cost[i][j] = min(cost[i][j],cost[i][k]+cost[k][j]);
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
    shortest_path();
    cin >> t;
    for (int i=1;i<=V;i++) cout << cost[t][i] << ' ';
    return 0;
}
#include <iostream>
#include <string>
#include <vector>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,f,V,E,t;
int d[MAX_N],cost[MAX_N][MAX_N];
bool used[MAX_N];
int pre[MAX_N];

void shortest_path(int s){ //O(V^2)
    fill(d,d+V+1,INF);
    fill(used,used+V+1,false);
    fill(pre, pre+V+1,-1);
    d[s] = 0;
    while (true){
        int v = -1;
        
        for (int u = 1;u <= V; u++){
            if (!used[u] && (v == -1 || d[u]<d[v])) v=u;
        }
        
        if (v == -1) break;
        used[v] = true;
        
        for (int u = 1; u <= V; u++){
            if (d[u]>d[v]+cost[v][u]){
                d[u] = d[v]+cost[v][u];
                pre[u] = v;
            }
        }
    }
}

void init(){
    cin >> V >> E;
    for (int i=1;i<=V;i++){
        for (int j=1;j<=V;j++){
            cost[i][j]=INF;
            cost[j][i]=INF;
        }
    }
    
    for (int i=0;i<E;i++){
        cin >> v >> u >> w;
        cost[v][u] = w;
        cost[u][v] = w;
    }
    
}

vector<int> get_path(int t){
    vector<int> path;
    for (; t!=-1; t = pre[t]) path.push_back(t);
    reverse(path.begin(),path.end());
    return path;
}

int main(){
    init();
    cout<<"DONE!\n";
    cin >> t;
    shortest_path(t);
    for (int i=1;i<=V;i++) cout << d[i] << ' ';
    cin >> t;
    vector<int> ans = get_path(t);
    for(int i=0; i<ans.size(); i++){
        cout << ans[i] << ' ';
    }
    return 0;
}
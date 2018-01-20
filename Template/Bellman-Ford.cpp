#include <iostream>
#include <string>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,f,V,E,t;
int d[MAX_N];

struct edge { int from,to,cost;};
edge es[MAX_N]; //edge set

void shortest_path(int s){ //O(VE)
    for (int i=1; i<=V; i++) d[i] = INF;
    d[s] = 0;
    while (true){
        bool update = false;
        for (int i=0; i<2*E; i++){
            edge e = es[i];
            if (d[e.from]!=INF && d[e.to]>d[e.from]+e.cost){
                d[e.to] = d[e.from] + e.cost;
                update = true;
            }
        }
        if (!update) break;
    }
}

bool find_negetive_loop(){
    memset(d,0,sizeof(d));
    for (int i=0;i<V;i++){
        for (int j=0;j<2*E;j++){
            edge e = es[j];
            if (d[e.to]>d[e.from]+e.cost){
                d[e.to] = d[e.from] + e.cost;
                if (i == V-1) return true;
            }
        }
    }
    return false;
}

void init(){
    cin >> V >> E;
    for (int i=0;i<E;i++){
        cin >> v >> u >> w;
        es[2*i].to = u;
        es[2*i].from = v;
        es[2*i].cost = w;
        es[2*i+1].to = v;
        es[2*i+1].from = u;
        es[2*i+1].cost = w;
    }
    
}

int main(){
    init();
    cout<<"DONE!\n";
    cin >> t;
    if (!find_negetive_loop()){
        shortest_path(t);
        for (int i=1;i<=V;i++) cout << d[i] << ' ';
    } else { // 一个负权值也是负环吗？
        cout << "Exist Negative Loop!";
    }
    return 0;
}
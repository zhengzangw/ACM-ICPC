#include <iostream>
#include <queue>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,V,E,t,nume,s;
int smin[MAX_N],cot[MAX_N];
bool h[MAX_N];

struct Edge
{
    int v,f,nxt;
};
Edge es[MAX_N*2+10];
int g[MAX_N + 10];

void addedge(int u, int v, int c){
    es[++nume].v=v;
    es[nume].f=c;
    es[nume].nxt=g[u];
    g[u]=nume;
}

void init(){
    cin >> V >> E;
    fill(g,g+V+1,0);
    nume=1;
    
    for (int i=0;i<E;i++){
        cin >> v >> u >> w;
        addedge(v,u,w);
        addedge(u,v,w);
    }
    
}

int shortest_path(int s){//O(ke) without SLF LLL
    fill(smin,smin+V+1,INF);
    fill(h,h+V+1,false);
    fill(cot,cot+V+1,0);
    smin[s] = 0;
    queue<int> que;
    while (!que.empty()) que.pop();
    que.push(s);
    h[s] = true; cot[s]++;
    
    while (!que.empty()){
        int u = que.front();
        que.pop();
        for (int v=g[u];v;v=es[v].nxt){
            if (smin[u]+es[v].f < smin[es[v].v]){
                smin[es[v].v] = smin[u] + es[v].f;
                if (!h[es[v].v]){
                    h[es[v].v]=true;
                    que.push(es[v].v);
                    cot[es[v].v]++;
                    if (cot[es[v].v]>V) return 1;
                }
            }
        }
        h[u] = false;
    }
    return 0;
}

int main(){
    init();
    cout<<"DONE!\n";
    cin >> s;
    if (shortest_path(s)) cout << "Existing Minus Loop!\n";
    else
        for (int i=1;i<=V;i++){
            cout << smin[i] << ' ';
        }
    
    return 0;
}
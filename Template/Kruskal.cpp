#include <iostream>
#include <string>
using namespace std;
const int MAX_N=1000, INF=1<<26;
int v,u,w,f,V,E,t;
int d[MAX_N];

struct edge { int from,to,cost;};
edge es[MAX_N]; //edge set

int par[MAX_N],ran[MAX_N];

void init_union_find(int n){
    for (int i = 0;i<n;i++){
        par[i] = i;
        ran[i] = 0;
    }
}

int find(int x){
    if (par[x] == x){
        return x;
    } else {
        return (par[x] = find(par[x]));
    }
}

int same(int x, int y){
    return find(x) == find(y);
}

void unite(int x, int y){
    x = find(x);
    y = find(y);
    if (x == y) return;
    
    if (ran[x]<ran[y]) par[x]=y;
    else {
        par[y]=x;
        if (ran[x]==ran[y]) ran[x]++;
    }
}

bool comp(const edge& e1,const edge& e2){
    return e1.cost<e2.cost;
}

int MST(){ //O(E*log(V))
    sort(es,es+E,comp);
    init_union_find(V);
    int res = 0;
    for (int i=0;i<E;i++){
        edge e = es[i];
        if (!same(e.from,e.to)){
            unite(e.from,e.to);
            res += e.cost;
        }
    }
    return res;
}

void init(){
    cin >> V >> E;
    for (int i=0;i<E;i++){
        cin >> v >> u >> w;
        es[i].to = u;
        es[i].from = v;
        es[i].cost = w;
    }
    
}

int main(){
    init();
    cout<<"DONE!\n";
    cout<<MST();
    return 0;
}
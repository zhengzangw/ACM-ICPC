#include <iostream>
using namespace std;
const int MAX_N = 1<<20;
int par[MAX_N],ran[MAX_N];
int n,m,t;
char c;

void init(int n){
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


int main(){
    cin >> n;
    init(100);
    for (int i=0;i<n;i++){
        cin >> c >> m >> t;
        if (c == 'q') {
            if (same(m,t)) cout << "Yes\n";
            else cout << "No\n";
        }
        if (c == 'u') {
            unite(m,t);
            cout << "Done!\n";
        }
    }
}
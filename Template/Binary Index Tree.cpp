#include <iostream>
using namespace std;
const int MAX_N = 1 << 26;
int n,m,t;
int bit[MAX_N],a[MAX_N];

int sum(int i){
    int s = 0;
    while (i>0){
        s += bit[i];
        i = i&(i-1); //or i-=i&-i
    }
    return s;
}

void add(int i, int x){
    while (i<=n){
        bit[i] += x;
        i+= i & -i;
    }
}

int main(){
    cin >> n;
    for (int i=1;i<=n;i++){
        cin >> t;
        add(i,t);
    }
    
    int s,t;
    char c;
    while (true){
        cin >> c;
        if (c=='e') break;
        if (c=='q'){
            cin >> s >> t;
            cout << sum(t)-sum(s-1) << '\n';
        }
        if (c=='s'){
                cin >> s >> t;
                add(s,t);
                cout << "Done with " << sum(s)-sum(s-1) << '\n';
        }
    }
    return 0;
}
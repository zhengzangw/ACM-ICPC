#include <iostream>
#include <set>
using namespace std;
int n,m,t;

int main(){
    cin >> n;
    
    set<int> s;
    //multiset<int> s can store repeated input
    
    for (int i=0;i<n;i++){
        cin >> t;
        s.insert(t);
    }
    
    set<int>::iterator ite;
    
    cin >> t;
    char c;
    while (t--){
        cin >> c >> m;
        if (c == 'd') s.erase(m);
        if (c == 'q') {
            // Way 1
            /* ite = s.find(m);
             if (ite == s.end()) cout << "Not Found!\n";
             else cout << "IN!\n"; */
            
            // Way 2
            if (s.count(m) != 0) cout << "Found!\n";
            else cout << "Not Found\n";
        }
        if (c == 'a') {
            for (ite = s.begin(); ite!=s.end();++ite){
                cout << *ite << ' ';
            }
            cout << '\n';
        }
    }
    return 0;
}
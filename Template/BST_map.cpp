#include <iostream>
#include <map>
#include <cstring>
using namespace std;
int n,k,t;
char s[100][100];

int main(){
    cin >> n;
    
    map<int , const char*> m;
    // <key, value>
    //multimap<int ,  const char*> m can store repeated input
    
    for (int i=0;i<n;i++){
        cin >> t >> s[i];
        m.insert(make_pair(t,s[i]));
    }
    
     map<int, const char*>::iterator ite;

    cin >> t;
    char c;
    while (t--){
        cin >> c >> k;
        if (c == 'd') m.erase(k);
        if (c == 'q') {
            // Way 1
            ite = m.find(k);
            if (ite == m.end()) cout << "Not Found!\n";
            else cout << ite->second << '\n';
            
            // Way 2
            /*if (m.count(k) != 0) cout << "Found!\n";
             else cout << "Not Found\n";*/
        }
        if (c == 'a') {
            for (ite = m.begin(); ite!=m.end();++ite){
                cout << ite->first << ' ' << ite->second << '\n';
            }
            cout << '\n';
        }
    }
    //revise: m[key] = new_value;
    return 0;
}
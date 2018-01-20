#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
using namespace std;

int kmp(string patten,string text){
    int n= patten.size();
    vector <int> next(n+1,0);
    for (int i=1;i<n;++i){
        int j=i;
        while (j>0){
            j = next[j];
            if (patten[i]==patten[j]) {
                next[i+1] = j+1;
                break;
            }
        }
    }
    
    int m = text.size();
    for (int i = 0, j = 0; i < m; ++i){
        if (j<n && text[i] == patten[j]) {j++;}
        else {
            while (j>0) {
                j = next[j];
                if (text[i] == patten[j]){
                    j++;
                    break;
                }
            }
        }
        
        if (j == n) return 1;
    }
    
    return 0;
    
}

int main(){
    cin.sync_with_stdio(false);
    int max_len,T,n;
    int flag;
    string tmp[200000];
    string max_s;
    cin >> T;
    while (T--){
        cin >> n;
        max_len = 0;
        for (int i=1;i<=n;i++){
            cin >> tmp[i];
            if (tmp[i].size() > max_len){
                max_len = tmp[i].size();
                max_s = tmp[i];
            }
        }
        flag = 1;
        for (int j=1; j<=n; j++){
            if (kmp(tmp[j],max_s) == 0){
                flag = 0;
                break;
                }
            }
        
        
        if (flag) cout << max_s << "\n";
        else cout << "No\n";
    }
}
#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N=1<<17;
int n,m,t,ans;
int smin[MAX_N][MAX_N];
int L,R;

int main(){
    cin >> n;
    for (int i=0;i<n;i++) cin >> smin[i][0];
    
    m = ceil(log2(n));
    for (int j=1; j<m;j++)
        for (int i=0; i+ (1<<j)-1 < n;i++)
            smin[i][j]= min(smin[i][j-1],smin[i+(1<<(j-1))][j-1]);
    
    cin >> L >> R;
    int k = log2(R-L+1);
    ans = min(smin[L][k],smin[R-(1 << k) + 1][k]);
    cout << ans;
    
    
    return 0;
}
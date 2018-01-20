#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N=1<<10;
int n,m,t,ans,is,ib;
int f[MAX_N],a[MAX_N],pres[MAX_N],preb[MAX_N];
int b[MAX_N],S[MAX_N],pret[MAX_N];

void init(){
    cin >> n;
    for (int i=0;i<n;i++) cin >> a[i];
}

int LSS(){
    for (int i=1;i<n;i++){
        f[i]=1;
        for (int j=0;j<i;j++)
            if (a[j]>=a[i] && f[j]+1>f[i]){
                f[i]=f[j]+1;
                pres[i]=j;
            }
    }
    
    int maxn=0;
    is=0;
    for (int i=0;i<n;i++) {
        if (maxn<f[i]) is=i;
        maxn=max(maxn,f[i]);
    }
    return maxn;
}

int BSe(int x,int y,int v){
    while(x<=y){
        int mid = x+(y-x)/2;
        if(S[mid]<=v) x=mid+1;
        else y=mid-1;
    }
    return x;
}

int LBS_fast(){
    int ans=0;
    for (int i=0;i<n;i++){
        int pos = BSe(1,i+1,a[i]);
        f[i] = pos;
        if (a[i]<S[pos]) {
            preb[i]=pret[pos-1];
            pret[pos]=i;
        }
        S[pos] = min(S[pos],a[i]);
        if (f[i]>ans) ib=i;
        ans = max(ans,f[i]);
    }
    return ans;
}

int LBS(){
    for (int i=1;i<n;i++){
        f[i]=1;
        for (int j=0;j<i;j++)
            if (a[j]<=a[i] && f[j]+1>f[i]){
                f[i]=f[j]+1;
                preb[i]=j;
            }
    }
    
    int maxn=0;
    for (int i=0;i<n;i++) {
        if (maxn<f[i]) ib=i;
        maxn=max(maxn,f[i]);
    }
    return maxn;
}

int main(){
    init();
    
    fill(f,f+n,0);
    fill(pres,pres+n,-1);
    f[0]=1;
    int s=LSS();
    
    fill(f,f+n,0);
    fill(preb,preb+n,-1);
    fill(pret,pret+n,-1);
    fill(b,b+n,0);
    fill(S,S+2*n,1<<20);
    int b=LBS_fast();
    
    ans = max(s,b);
    int temp[MAX_N];
    cout <<"The length of longest monotone sequence is "<< ans << '\n';
    if (s>=b) {
        cout << "A not increased one:\n";
        int k=0;
        for (int i=is;i>-1;i=pres[i]) temp[k++]=a[i];
        for (int i=k-1;i>=0;i--) cout << temp[i] << ' ';
        cout << '\n';
    }
    if (s<=b) {
        cout << "A not decreased one:\n";
        int k=0;
        for (int i=ib;i>-1;i=preb[i]) temp[k++]=a[i];
        for (int i=k-1;i>=0;i--) cout << temp[i] << ' ';
        cout << '\n';
    }
    return 0;
}

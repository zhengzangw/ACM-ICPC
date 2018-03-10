#include <iostream>
#include <cstdio>
#include <cmath>
#include <vector>
using namespace std;
const int MAX_N = 1000;
int n,m,t,B;
int a[MAX_N];
int I[MAX_N],J[MAX_N],K[MAX_N],nums[MAX_N];//查询i到j升序排列第k个数

void solve(){
    vector<int> bucket[n/B+1];
    for (int i=0;i<n;i++){
        bucket[i/B].push_back(a[i]);
        nums[i] = a[i];
    }
    sort(nums,nums+n);
    for (int i=0;i<n/B;i++) sort(bucket[i].begin(),bucket[i].end());
    //最后多余的桶没有排序
    for (int i=0;i<m;i++){
        int l=I[i], r=J[i]+1, k=K[i];//l正好，r指下一个
        int lb=-1,ub=n-1;
        while (ub-lb>1){
            int md = (lb+ub)/2;
            int x = nums[md];
            int tl = l, tr = r, c = 0;
            
            while (tl<tr && tl % B != 0) if (a[tl++]<=x) c++;
            while (tl<tr && tr % B != 0) if (a[--tr]<=x) c++;
            
            while (tl<tr){
                int b = tl/B;
                c += upper_bound(bucket[b].begin(),bucket[b].end(),x) - bucket[b].begin();
                tl += B;
            }
            
            if (c>=k) ub = md;
            else lb=md;
        }
        
        printf("%d\n",nums[ub]);
    }
}
int main(){
    cin >> n >> m;
    for (int i=0;i<n;i++) cin >> a[i];
    for (int i=0;i<m;i++) {
        cin >> I[i] >> J[i] >> K[i];
        I[i]--; J[i]--;
    }
    B = floor(sqrt(n));
    solve();
    return 0;
}

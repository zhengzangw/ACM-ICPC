#include <iostream>
#include <cmath>
using namespace std;
const int MAX_N=1<<17;
int n,m,t;
struct segtree{
    int ma,mi,sum,lazy;
};
segtree seg[MAX_N],a[MAX_N];

void PushUp(int node)
{
    seg[node].mi = min(seg[2*node+1].mi,seg[2*node].mi);
    seg[node].ma = max(seg[2*node+1].ma,seg[2*node].ma);
    seg[node].sum = seg[2*node+1].sum+seg[2*node].sum;
}

void build(int node,int begin, int end)
{
    if (begin == end) seg[node] = a[begin];
    else{
        build(2*node,begin,(begin+end)/2);
        build(2*node+1,(begin+end)/2+1,end);
        //不重叠线段树
        
        PushUp(node);
    }
}

void PushDown(int node,int left,int right,int m){
    if (seg[node].lazy){
        seg[node<<1].lazy+=seg[node].lazy;
        seg[(node<<1)+1].lazy+=seg[node].lazy;
        seg[node<<1].sum += seg[node].lazy*(m-left+1);
        seg[(node<<1)+1].sum += seg[node].lazy*(right-m);
        seg[node<<1].mi += seg[node].lazy*(m-left+1);
        seg[(node<<1)+1].mi += seg[node].lazy*(right-m);
        seg[node<<1].ma += seg[node].lazy*(m-left+1);
        seg[(node<<1)+1].ma += seg[node].lazy*(right-m);
        seg[node].lazy=0;
    }
}

int query_min(int node, int begin, int end, int left, int right)
{
    if (begin >= left && end <= right) return seg[node].mi;
    PushDown(node,begin,end,m);
    int m=(begin+end) >> 1;
    
    int p1=1<<26,p2=1<<26;
    if (left<=m) p1=query_min(2*node,begin,m,left,right);
    if (m<right) p2=query_min(2*node+1,m+1,end,left,right);
    return min(p1,p2);
}

int query_max(int node, int begin, int end, int left, int right)
{
    if (begin >= left && end <= right) return seg[node].ma;
    PushDown(node,begin,end,m);
    int m=(begin+end) >> 1;
    
    int p1=-1,p2=-1;
    if (left<=m) p1=query_max(2*node,begin,m,left,right);
    if (m<right) p2=query_max(2*node+1,m+1,end,left,right);
    return max(p1,p2);
}

int query_sum(int node, int begin, int end, int left, int right)
{
    if (begin >= left && end <= right) return seg[node].sum;
    PushDown(node,begin,end,m);
    int m=(begin+end) >> 1;
    
    int p1=0,p2=0;
    if (left<=m) p1=query_sum(2*node,begin,m,left,right);
    if (m<right) p2=query_sum(2*node+1,m+1,end,left,right);
    return p1+p2;
}

void update_node (int node, int begin, int end, int ind, int add)
{
    if (begin == end)
    {
        seg[node].mi = add;
        seg[node].ma = add;
        seg[node].sum = add;
        return;
    }

    int m = (begin+end) >> 1;
    PushDown(node,begin,end,m);
    if (ind <= m) update_node(node*2,begin,m,ind,add);
    else update_node(node*2+1, m+1, end, ind, add);
    PushUp(node);
}

void change(int node, int begin, int end, int left, int right, int add)
{
    if (left <= begin && end <= right) {
        seg[node].lazy += add;
        seg[node].sum += add*(end-begin+1);
        seg[node].mi += add;
        seg[node].ma += add;
        return;
    }
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    if (left<=m) change(2*node,begin,m,left,right,add);
    if (m<right) change(2*node+1,m+1,end,left,right,add);
    PushUp(node);
}

int main(){
    cin >> n;
    for (int i=1;i<=n;i++){
        cin >> t;
        a[i].ma=t;
        a[i].mi=t;
        a[i].sum=t;
    }
    build(1,1,n);
    
    int s,t,w;
    char c;
    while (true){
        cin >> c;
        if (c=='e') break;
        if (c=='i'){
            cin >> s >> t;
            cout << query_min(1,1,n,s,t) << '\n';
        }
        if (c=='a'){
			cin >> s >> t;
            cout << query_max(1,1,n,s,t) << '\n';
        }
        if (c=='q'){
        	cin >> s >> t;
            cout << query_sum(1,1,n,s,t) << '\n';
        }
        if (c=='s'){
            cin >> s >> t;
            update_node(1,1,n,s,t);
            cout << "Done with " << query_min(1,1,n,s,s) << '\n';
        }
        if (c=='c'){
            cin >> s >> t >> w;
            change(1,1,n,s,t,w);
            cout << "Done!\n";
        }
    }
    return 0;
}
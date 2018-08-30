#ifndef Seg_Tree_H
#define Seg_Tree_H
#include <bits/stdc++.h>
using namespace std;
const int MAXN=1e5;
struct Segnode{
    int ma,mi,sum,add_lazy,mul_lazy;
    Segnode(){}
    Segnode(int n):ma(n),mi(n),sum(n),add_lazy(0){}
};
Segnode seg[MAXN];

void PushUp(int node)
{
    seg[node].mi = min(seg[2*node+1].mi,seg[2*node].mi);
    seg[node].ma = max(seg[2*node+1].ma,seg[2*node].ma);
    seg[node].sum = seg[2*node+1].sum+seg[2*node].sum;
}

void PushDown(int node,int left,int right,int m){
    if (seg[node].add_lazy){ //||seg[node].mul!=1
        int lson = node<<1, rson = (node<<1)+1;
        /*
        seg[lson].mull = seg[lson].mull * seg[node].mull % mod;
        seg[rson].mull = seg[rson].mull * seg[node].mull %mod;
        seg[lson].lazy = seg[lson].lazy * seg[node].mull %mod;
        seg[rson].lazy = seg[rson].lazy * seg[node].mull %mod;
        seg[lson].sum  = seg[lson].sum  * seg[node].mull %mod;
        seg[rson].sum  = seg[rson].sum  * seg[node].mull %mod;
        seg[node].mull=1;
        */
        seg[lson].add_lazy+=seg[node].add_lazy;
        seg[rson].add_lazy+=seg[node].add_lazy;
        seg[lson].sum += seg[node].add_lazy*(m-left+1);
        seg[rson].sum += seg[node].add_lazy*(right-m);
        seg[lson].mi += seg[node].add_lazy*(m-left+1);
        seg[rson].mi += seg[node].add_lazy*(right-m);
        seg[lson].ma += seg[node].add_lazy*(m-left+1);
        seg[rson].ma += seg[node].add_lazy*(right-m);
        seg[node].add_lazy=0;
    }
}

//初始化
void build(int node,int begin, int end)
{
    if (begin == end) seg[node] = Segnode(0);
    else{
        build(2*node,begin,(begin+end)/2);
        build(2*node+1,(begin+end)/2+1,end);
        //不重叠线段树
        PushUp(node);
    }
}

//查询
int query_min(int node, int begin, int end, int left, int right)
{
    //if (left>right) return 0;
    if (begin >= left && end <= right) return seg[node].mi;
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    
    int p1=INT_MAX,p2=INT_MAX;
    if (left<=m) p1=query_min(2*node,begin,m,left,right);
    if (m<right) p2=query_min(2*node+1,m+1,end,left,right);
    return min(p1,p2);
}

int query_max(int node, int begin, int end, int left, int right)
{
    //if (left>right) return 0;
    if (begin >= left && end <= right) return seg[node].ma;
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    
    int p1=-INT_MAX,p2=-INT_MAX;
    if (left<=m) p1=query_max(2*node,begin,m,left,right);
    if (m<right) p2=query_max(2*node+1,m+1,end,left,right);
    return max(p1,p2);
}

int query_sum(int node, int begin, int end, int left, int right)
{
    //if (left>right) return 0;
    if (begin >= left && end <= right) return seg[node].sum;
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    
    int p1=0,p2=0;
    if (left<=m) p1=query_sum(2*node,begin,m,left,right);
    if (m<right) p2=query_sum(2*node+1,m+1,end,left,right);
    return p1+p2;
}

//单点更新
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

//区间更新
void interval_add(int node, int begin, int end, int left, int right, int add)
{
    if (left <= begin && end <= right) {
        seg[node].add_lazy += add;
        seg[node].sum += add*(end-begin+1);
        seg[node].mi += add;
        seg[node].ma += add;
        return;
    }
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    if (left<=m) interval_add(2*node,begin,m,left,right,add);
    if (m<right) interval_add(2*node+1,m+1,end,left,right,add);
    PushUp(node);
}
void interval_mul(int node, int begin, int end, int left, int right, int mul){
    if (left<=begin && end <= right){
        seg[node].mul_lazy = (seg[node].mul_lazy*mul);
        seg[node].add_lazy = (seg[node].add_lazy*mul);
        seg[node].sum  = (seg[node].sum *mul);
        return;
    }
    int m=(begin+end) >> 1;
    PushDown(node,begin,end,m);
    if (left<=m) mul(2*node,begin,m,left,right,mul);
    if (m<right) mul(2*node+1,m+1,end,left,right,mul);
    PushUp(node);
}

int discrete(int n,int a[]){
    int d[MAXN];
    memcpy(d,a,sizeof(d));
    sort(d,d+n);
    int len = int(unique(d,d+n) - d);
    for (int i=0;i<n;++i) a[i] =int(lower_bound(d,d+n,a[i]) - d + 1);
    return len;
}
#endif
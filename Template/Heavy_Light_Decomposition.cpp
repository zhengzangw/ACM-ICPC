#include <cstdio>
#include <algorithm>
#include <cstring>
using namespace std;
const int MAX_N = 20000;
int n,m,t;
int root,edge,a,b,c,z;
int tree[MAX_N],d[MAX_N][3],first[MAX_N],dep[MAX_N],w[MAX_N],
top[MAX_N],son[MAX_N],siz[MAX_N],fa[MAX_N],lazy[MAX_N];
char ch[10];
struct Tedge {int b,next;} e[MAX_N*2];

void insert (int a,int b,int c)
{
    e[++edge].b=b;
    e[edge].next = first[a];
    first[a]=edge;
}

void dfs(int v)//calculate fa,son,dep,siz
{
    siz[v]=1; son[v]=0;
    for (int i=first[v];i>0;i=e[i].next)
        if (e[i].b!=fa[v])
        {
            fa[e[i].b]=v;
            dep[e[i].b]=dep[v]+1;
            dfs(e[i].b);
            if (siz[e[i].b]>siz[son[v]]) son[v]=e[i].b;
            siz[v] += siz[e[i].b];
        }
}


void PushDown(int root,int left,int right,int m){
    if (lazy[root]){
        lazy[root*2] += lazy[root];
        lazy[root*2+1] += lazy[root];
        tree[root*2] += lazy[root];
        tree[root*2+1] += lazy[root];
        lazy[root] = 0;
    }
}

void build_tree(int v, int tp)//calculate w,top
{
    w[v] = ++z; top[v] = tp;
    if (son[v]!=0) build_tree(son[v],top[v]);
    for (int i=first[v];i>0;i=e[i].next)
        if (e[i].b!=son[v] && e[i].b!= fa[v])
            build_tree(e[i].b,e[i].b);
}

void update(int root, int left, int right, int ind,int add)//segment tree
{
    if (ind > right || left > ind) return;
    if (left==right)
    { tree[root]=add; return;}
    int mid = (left+right)/2, ls = root * 2, rs = ls+1;
    PushDown(root,left,right,mid);
    update(ls,left,mid,ind,add);
    update(rs,mid+1,right,ind,add);
    tree[root] = max(tree[ls],tree[rs]);
}

void add_all(int root, int left, int right, int L, int R, int add)//segment tree
{
    if (L>right || R<left) return;
    if (L<=left && right<=R) {
        tree[root] += add;
        lazy[root] += add;
        return;
    }
    int mid = (left+right)/2, ls = root*2, rs=ls+1;
    PushDown(root,left,right,mid);
    add_all(ls,left,mid,L,R,add);
    add_all(rs,mid+1,right,L,R,add);
}

void update_all(int va, int vb, int add)
{
    int f1=top[va],f2=top[vb];
    while (f1!=f2)
    {
        if (dep[f1]<dep[f2])
        {swap(f1,f2); swap(va,vb);}
        add_all(1,1,z,w[f1],w[va],add);
        va=fa[f1]; f1=top[va];
    }
    if (va == vb) return;
    if (dep[va]>dep[vb]) swap(va,vb);
    add_all(1,1,z,w[son[va]],w[vb],add);
}

int maxi(int root, int left,int right, int L, int R)//segment tree
{
    if (L>right || R<left) return 0;
    if (L<=left && right<=R) return tree[root];
    int mid = (left+right)/2, ls = root*2, rs=ls+1;
    PushDown(root,left,right,mid);
    return max(maxi(ls,left,mid,L,R),maxi(rs,mid+1,right,L,R));
}

int find(int va, int vb)
{
    int f1=top[va],f2=top[vb],tmp=0;
    while (f1!=f2)
    {
        if (dep[f1]<dep[f2])
        {swap(f1,f2); swap(va,vb);}
        tmp = max(tmp,maxi(1,1,z,w[f1],w[va]));
        va=fa[f1]; f1=top[va];
    }
    if (va == vb) return tmp;
    if (dep[va]>dep[vb]) swap(va,vb);
    return max(tmp,maxi(1,1,z,w[son[va]],w[vb]));
}

void init()
{
    scanf("%d",&n);
    root = (n+1)/2;
    fa[root] = z = dep[root] = edge =0;
    memset(siz,0,sizeof(siz));
    memset(first,0,sizeof(first));
    memset(tree,0,sizeof(tree));
    memset(lazy,0,sizeof(lazy));
    for (int i=1;i<n;i++){
        scanf("%d%d%d",&a,&b,&c);
        d[i][0]=a; d[i][1]=b; d[i][2]=c;
        insert(a,b,c);
        insert(b,a,c);
    }
    
    dfs(root);
    build_tree(root,root);
    for (int i=1;i<n;i++){
        if (dep[d[i][0]] > dep[d[i][1]]) swap(d[i][0],d[i][1]);
        update(1,1,z,w[d[i][1]],d[i][2]);
    }
    
}

inline void read()
{
    ch[0]=' ';
    while (ch[0]<'A'||ch[0]>'Q') scanf("%s",&ch);
}

void solve()
{
    for (read();ch[0]!='D';read()){
        scanf("%d%d",&a,&b);
        if (ch[0]=='Q') printf("%d\n",find(a,b));
        if (ch[0]=='C') update(1,1,z,w[d[a][1]],b);
        if (ch[0]=='A') {
            scanf("%d",&c);
            update_all(a,b,c);
        }
    }
}

int main(){
    scanf("%d",&t);
    while (t--)
    {
        init();
        solve();
    }
    return 0;
}

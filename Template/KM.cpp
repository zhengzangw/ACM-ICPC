#include<cstdio>
#include<cstring>
#include<queue>
#include<cmath>

using namespace std;
const int Maxn = 1100;
const int INF = 1 << 26;
int n,m,t;
int w[Maxn][Maxn],linky[Maxn],lx[Maxn],ly[Maxn],visx[Maxn],visy[Maxn];
int lack;

bool find(int x)
{
    visx[x] = true;
    for (int y=1;y<=n;y++){
        if (visy[y]) continue;
        int t = lx[x]+ly[y]-w[x][y];
        if (!t) {
            visy[y] = true;
            if (!~linky[y]||find(linky[y])){
                linky[y] = x;
                return true;
            }
        } else lack = min(lack,t);
    }
    return false;
}

int KM(){//O(n3)
    memset(linky,-1,sizeof(linky));
    for (int i=1;i<=n;++i)
        for (int j=1;j<=m;++j)
            lx[i] = max(lx[i],w[i][j]);
    
    for (int x=1;x<=n;++x)
        while (true){
            memset(visx,0,sizeof(visx));
            memset(visy,0,sizeof(visy));
            lack = INF;
            if (find(x)) break;
            for (int i=1;i<=n;++i){
                if (visx[i]) lx[i]-=lack;
                if (visy[i]) ly[i]+=lack;
            }
        }
    
    int ans=0;
    for (int i=1;i<=n;i++) ans+=lx[i];
    for (int i=1;i<=m;i++) ans+=ly[i];
    return ans;
}

int main(){
    scanf("%d%d",&n,&m);
    for (int i=1;i<=n;i++)
        for (int j=1;j<=m;j++)
            scanf("%d",&w[i][j]);
    /*for (int i=1;i<=n;i++)
     for (int j=1;j<=m;j++)
     w[i][j]=-w[i][j];*/
    //min printf -KM()
    printf("%d\n",KM());
    return 0;
}

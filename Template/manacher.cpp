#include <iostream>
#include <cstdio>
#include <algorithm>
using namespace std;
const int maxn=100;
int n,i,j;
char s[maxn],tmp[maxn];
int L[maxn];

int init(char *st)
{
    int len=(int) strlen(s);
    tmp[0]='@';
    for (int i=1;i<=2*len;i+=2){
        tmp[i]='#';
        tmp[i+1]=s[i/2];
    }
    tmp[2*len+1]='#';
    tmp[2*len+1]='$';
    tmp[2*len+3]=0;
    return 2*len+1;
}

int manacher(char *st)
{
    int len =(int) strlen(tmp)-2;
    int mx=0,ans=0,po=0;
    for (int i=1;i<=len;i++)
    {
        if (mx>i) L[i]=min(L[2*po-i],mx-i);
        else L[i]=1;
        
        while (tmp[i-L[i]]==tmp[i+L[i]]) L[i]++;
        
        if (L[i]+i>mx)
        {
            mx = L[i]+1;
            po = i;
        }
        ans = max(ans,L[i]);
    }
    return ans-1;
}

int main(){
    scanf("%s",&s);
    init(s);
    cout << manacher(tmp);
    return 0;
}

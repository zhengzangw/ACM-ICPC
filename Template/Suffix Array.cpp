#include <iostream>
#include <cstdio>
#include <queue>
#include <algorithm>
using namespace std;
#define MAX 26
int n,m,cnt,k,h;
int sa[100],tmp[100],rnk[100],lcp[100];

bool compare_sa(int i, int j)
{
    if (rnk[i]!=rnk[j]) return rnk[i]<rnk[j];
    
    int ri = i+k<=n ? rnk[i+k] : -1;
    int rj = j+k<=n ? rnk[j+k] : -1;
    return ri < rj;
}

void Manber_Myers(string s) //O(nlog2n)
{
    n = s.length();
    for (int i=0;i<=n;i++)
    {
        sa[i] = i;
        rnk[i] = i < n ? s[i] : -1;
    }
    
    for (k=1;k<=n;k*=2)
    {
        sort(sa,sa+n+1,compare_sa);
        tmp[sa[0]] = 0;
        for (int i = 1;i<=n;i++)
            tmp[sa[i]] = tmp[sa[i-1]] + (compare_sa(sa[i-1],sa[i]) ? 1:0);
        for (int i = 0; i<=n; i++) rnk[i] = tmp[i];
    }
}

int contain(string s,string t)
{
    int a=0,b=s.length();
    while (b-a>1){
        int c=(a+b)/2;
        if (s.compare(sa[c],t.length(),t)<0) a = c;
        else b = c;
    }
    if (s.compare(sa[b],t.length(),t)==0) return sa[b];
    else return -1;
}

void construct_lcp(string s)
{
    int n=s.length();
    for (int i=0;i<=n;i++) rnk[sa[i]]=i;

    int h = 0;
	lcp[0] = 0;
	for (int i=0;i<n;i++){
		int j = sa[rank[i]-1];
		if (h>0) h--;
		for (;j+h<n && i+h<n;h++){
			if (S[j+h]!=S[i+h]) break;
		}
		lcp[rank[i]-1] = h;
	}
}

int Longest_Common_String(string s, string t)
{
	int n=s.length();
	string s+='$'+t;

	Manber_Myers(s);
	construct_lcp(s);

	int ans = 0;
	for (int i = 0;i<s.length();i++)
	 ans = ((sa[i]<n)!=(sa[i]>n)) ? max(ans,lcp[i]);
	return ans;
}
//With rmq, we can calculate the longest perfix of any two of the suffix.

int main(int argc, char *argv[])
{
    string s,t;
    cin >> s >> t;
    Manber_Myers(s);
    cout << contain(s,t)+1 << '\n';
    construct_lcp(s);
    cin >> t;
    Longest_Common_String(s,t);
    return 0;
}


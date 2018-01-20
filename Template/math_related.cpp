#include <cstdio>
#include <vector>
typedef long long LL;
using namespace std;

void swap(LL a, LL b){
  LL t = b;
  b = a; a = t;
}

LL gcd(LL a,LL b){
  return b == 0 ? a:gcd(b,a % b);
}

LL extend_gcd( LL a, LL b, LL& x, LL& y){
  if (b==0) {
    x=1; y=0;
    return a;
  } else {
    LL r = extend_gcd(b,a % b,y,x);
    y -= x*(a/b);
    return r;
  }
}

LL ni(LL x, LL mod){ // only for (x,mod)=1
   LL xx,yy;
   if (gcd(x,mod)>1) return -1;
   extend_gcd(mod,x,xx,yy);
   return yy;
}

LL frac_mod(LL a,LL b,LL m){
  if (gcd(b,m)>1) return -1;
  return (ni(b,m) * a %m);
}

vector<LL> line_mod_equation(LL a, LL b, LL n){
  LL x,y;
  LL d = extend_gcd(a,b,x,y);
  vector <LL> ans;
  ans.clear();
  if (b%d == 0) {
    x %= n; x += n; x %= n;
    ans.push_back(x*(b/d)%n);
    for (LL i=1; i<d; ++i)
      ans.push_back((ans[0]+i*n/d)%n);
  }
  return ans;
}

LL pow(LL a, LL i, LL n){
  if (i==0) return 1;
  LL temp = pow(a,i>>1,n);
  temp = temp*temp % n;
  if (i & 1) temp = temp * a % n;
  return temp;
}

int main(){
  LL a,b,x,y,d,ni_x,mod,n;
  scanf("%lld%lld",&a,&b);
  if (a<b) swap(a,b);
  printf("%lld\n",gcd(a,b));
  d=extend_gcd(a,b,x,y);
  printf("%lldA+%lldB=%lld",x,y,d);
  
  scanf("%lld%lld",&ni_x,&mod);
  printf("%lld*ni_x=%lld",ni(ni_x,mod),mod);

  scanf("%lld%lld%lld",&a,&b,&mod);
  printf("%lld",frac_mod(a,b,mod));

  scanf("%lld%lld%lld",&a,&b,&n);
  vector <LL> ans = line_mod_equation(a,b,n);
  int i;
  for (i=0;i<ans.size();i++)
    printf("%lld\n",ans[i]);

  scanf("%lld%lld%lld",&a,&n,&mod);
  printf("%lld",pow(a,n,mod));

  return 0;
}
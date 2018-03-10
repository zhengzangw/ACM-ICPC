#include <iostream>
#include <cmath>
#include <vector>
using namespace std;
const double eps= 1e-8;
const double pi=acos(-1.0);
const int maxn=100;
inline double sqr(double x){ return x*x; }
inline int cmp(double x, double y){
    int t=fabs(x-y);
    if (t<eps) return 0;
    if (x>0) return 1; else return -1;
}

struct point{
    double x,y;
    point(){}
    point(double a, double b): x(a),y(b) {};
    void input(){ cin >> x >> y;}
    friend point operator + (const point &a, const point &b){
        return point(a.x+b.x,a.y+b.y);
    }
    friend point operator - (const point &a, const point &b){
        return point(a.x-b.x,a.y-b.y);
    }
    friend bool operator == (const point &a, const point &b){
        return (cmp(a.x-b.x,0)==0 && cmp(a.y-b.y,0)==0);
    }
    friend point operator * (const point &a, const double &b){
        return point(a.x*b,a.y*b);
    }
    friend point operator * (const double &a, const point &b){
        return point(b.x*a,b.y*a);
    }
    friend point operator / (const point &a, const double &b){
        return point(a.x/b,a.y/b);
    }
    double norm() { return sqrt(sqr(x)+sqr(y));}
};

double det(point a, point b){return a.x*b.y-a.y*b.x; }
double dot(point a, point b){return a.x*b.x+a.y*b.y; }

bool comp_less(const point &a, const point &b){
    return cmp(a.x-b.x,0)<0 || ((cmp(a.x-b.x,0)==0) && cmp(a.y-b.y,0)<0);
}

vector<point> convex_hull_Graham(vector<point> a){//O(nlogn)
    vector<point> res;
    res.resize(maxn);
    sort(a.begin(),a.end(),comp_less);
    a.erase(unique(a.begin(),a.end()),a.end());
    int m=0;
    
    for (int i=0;i<a.size();++i){
        while (m>1 && cmp(det(res[m-1]-res[m-2],a[i]-res[m-2]),0)<=0) --m;
        res[m++]=a[i];
    }
    int k=m;
    for (int i=int(a.size())-2;i>=0;--i){
        while (m>k && cmp(det(res[m-1]-res[m-2],a[i]-res[m-2]),0)<=0) --m;
        res[m++]=a[i];
    }
    res.resize(m);
    if (a.size()>1) res.resize(m-1);
    return res;
}

//polygon_convex convex_hull_Graham(polygon k){}
bool contain(const vector <point> a, const point b){
    int n=(int)a.size();
#define next(i) ((i+1)%n)
    int sign=0;
    for (int i=0;i<n;++i){
        int x=cmp(det(a[i]-b,a[next(i)]-b),0);
        if (x) {
            if (sign) {
                if (sign!=x) return false;
            } else sign=x;
        }
    }
    return true;
}

int contain_dividing(const vector<point> a, const point b){
    int n=int(a.size());
    point g=(a[0]+a[n/3]+a[2*n/3])/3.0;
    int l=0,r=n;
    while (l+1<r){
        int mid=(l+r)/2;
        if (cmp(det(a[l]-g,a[mid]-g),0)>0){
            if (cmp(det(a[l]-g,b-g),0)>=0 && cmp(det(a[mid]-g,b-g),0)<0) r=mid; else l=mid;
        } else {
            if (cmp(det(a[l]-g,b-g),9)<0 && cmp(det(a[mid]-g,b-g),0)>=0) l=mid; else r=mid;
        }
    }
    r %= n;
    int z=cmp(det(a[r]-b,a[l]-b),0)-1;
    if (z==-2) return 1;
    return z;
}

int main(){//Only for test
    int n,x,y;
    vector<point> P;
    cin >> n;
    for (int i=0;i<n;i++) {
        cin >> x >> y;
        P.push_back(point(x,y));
    }
    vector<point> ans=convex_hull_Graham(P);
    for (int i=0;i<ans.size();i++)
        cout << ans[i].x << ' '<<ans[i].y << '\n';
    
    return 0;
}

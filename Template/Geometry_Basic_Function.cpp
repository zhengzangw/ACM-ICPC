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

inline int gcd(int a, int b){
    if (a<b) swap(a,b);
    return b==0 ? a:gcd(b,a%b);
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

double dist(point a, point b){return (a-b).norm();}

point rotate_point(point p, double A){//逆时针
    double tx=p.x,ty=p.y;
    return point(tx*cos(A)-ty*sin(A),tx*sin(A)+ty*cos(A));
}

struct line{
    point a,b;
    bool segment;
    line() {}
    line(point x,point y): a(x),b(y),segment(1){}
    line(point x, point y, bool flag): a(x),b(y),segment(flag){}
    point line2point(){ return point(b-a);}
};

double dis_point_line(const point p,const line L){//点到直线距离
    point s=L.a,t=L.b;
    if (L.segment==0) return fabs(det(s-p,t-p)/dist(s,t));
    if (cmp(dot(p-s,t-s),0)==-1) return (p-s).norm();
    if (cmp(dot(p-t,s-t),0)==-1) return (p-t).norm();
    return fabs(det(s-p,t-p)/dist(s,t));//way 2 dist(point_proj_line(p,L),p);
}

point point_proj_line(const point p,const line L){//垂足
    point s=L.a,t=L.b;
    double r=dot(t-s,p-s)/dot(t-s,t-s);
    return s+r*(t-s);
}
bool parallel(line a,line b){//平行 including overlap
    return !cmp(det(a.a-a.b,b.a-b.b),0);
}

bool overlap(line a,line b){
    return parallel(a,b) && (cmp(det(a.a-b.a,a.a-b.b),0)==0);
}

bool pointinrec(point p, line L){ //点在矩形内
    if (L.segment==0) return true;
    if (min(L.a.x,L.b.x) <= p.x <= max(L.a.x,L.b.x) && min(L.a.y,L.b.y) <= p.y <= max(L.a.y,L.b.y))
        return true;
    return false;
}

bool point_on_line(point p, line L){//点在线上
    if (cmp(det(p-L.a,L.a-L.b),0)!=0) return false;
    if (L.segment == 0) return true;
    else if (cmp(dot(p-L.a,p-L.b),0)!=1) return true;
    else return false;
    // Way 2: return pointinrec(p,L);
}

bool recinrec(line l1, line l2){//矩形在矩形中
    if (l1.segment == 0 || l2.segment == 0) return true;
    point s,t;
    s.x=max(l1.a.x,l2.a.x); s.y=max(l1.a.y,l2.a.y);
    t.x=min(l1.b.x,l2.b.x); t.y=min(l1.b.y,l2.b.y);
    line L_rec(s,t,1);
    if (L_rec.a.x > L_rec.b.x || L_rec.a.y > L_rec.b.y) return false;
    return true;
}

bool line_make_point(line a, line b, point &res) {//交点判断
    if (overlap(a,b)) {res.x=1 << 20;res.y=1 << 20; return true;};
    if (parallel(a,b)) return false;
    //if (!recinrec(a,b)) return false;
    double s1 = det(a.a-b.a,b.b-b.a);
    double s2 = det(a.b-b.a,b.b-b.a);
    res = (s1*a.b-s2*a.a)/(s1-s2);
    if (point_on_line(res,a) && point_on_line(res,b)) return true;
    else return false;
    /* Way two 跨立测试
     if (a.segment == 1 && b.segment == 1){
     if (!recinrec(l1,l2)) return false;
     if (cross(l1.a-l2.a,l2.b-l2.a)*cross(l1.b-l2.a,l2.b-l2.a)>0) return false;
     if (cross(l2.a-l1.a,l1.b-l1.a)*cross(l2.b-l1.a,l1.b-l1.a)>0) return false;
     return true;
     }*/
}

line move_d(line a, double len){//法向移动
    point d=a.line2point();
    d = d/d.norm();
    d = rotate_point(d,pi/2);
    return line(a.a+d*len,a.b+d*len);
}

struct polygon {
    int n;
    point a[maxn];
    polygon(){}
    double perimeter(){
        double sum=0;
        a[n]=a[0];
        for (int i=0;i<n;i++)
            sum+=(a[i+1]-a[i]).norm();
        return sum;
    }
    double area() {
        double sum=0;
        a[n]=a[0];
        for (int i=0;i<n;i++) sum+=det(a[i+1],a[i]);
        return fabs(sum/2.);
    }
    int Point_In(point t){
        int num=0,i,d1,d2,k;
        a[n]=a[0];
        for (i=0;i<n;i++){
            line L(a[i],a[i+1],1);
            if (point_on_line(t,L)) return 2;
            k = cmp(det(a[i+1]-a[i],t-a[i]),0);
            d1 = cmp(a[i].y-t.y,0);
            d2 = cmp(a[i].x-t.x,0);
            if (k>0 && d1<=0 && d2>0) num++;
            if (k<0 && d2<=0 && d1>0) num--;
                }
        return num!=0;
    }
    point MassCenter(){//要求多边形点逆时针输入
        point ans=point(0,0);
        if (cmp(area(),0)==0) return ans;
        a[n]=a[0];
        for (int i=0;i<n;i++) ans = ans+(a[i]+a[i+1])*det(a[i+1],a[i]);
        return ans/area()/6.; // /2 /3
    }
    int Border_Point(){
        int num=0;
        a[n]=a[0];
        for (int i=0;i<n;i++)
            num+=gcd(abs(int(a[i+1].x-a[i].x)),abs(int(a[i+1].y-a[i].y)));
        return num;
    }
};

struct polygon_convex{
    vector <point> P;
    polygon_convex(int size=0){
        P.resize(size);
    }
};

bool comp_less(const point &a, const point &b){
    return cmp(a.x-b.x,0)<0 || ((cmp(a.x-b.x,0)==0) && cmp(a.y-b.y,0)<0);
}

polygon_convex convex_hull_Graham(vector<point> a){//O(nlogn)
    polygon_convex res(2*(int)a.size()+5);
    sort(a.begin(),a.end(),comp_less);
    a.erase(unique(a.begin(),a.end()),a.end());
    int m=0;
    
    for (int i=0;i<a.size();++i){
        while (m>1 && cmp(det(res.P[m-1]-res.P[m-2],a[i]-res.P[m-2]),0)<=0) --m;
        res.P[m++]=a[i];
    }
    int k=m;
    for (int i=int(a.size())-2;i>=0;--i){
        while (m>k && cmp(det(res.P[m-1]-res.P[m-2],a[i]-res.P[m-2]),0)<=0) --m;
        res.P[m++]=a[i];
    }
    res.P.resize(m);
    if (a.size()>1) res.P.resize(m-1);
    return res;
}

//polygon_convex convex_hull_Graham(polygon k){}
bool contain(const polygon_convex &a, const point &b){
    int n=(int)a.P.size();
#define next(i) ((i+1)%n)
    int sign=0;
    for (int i=0;i<n;++i){
        int x=cmp(det(a.P[i]-b,a.P[next(i)]-b),0);
        if (x) {
            if (sign) {
                if (sign!=x) return false;
            } else sign=x;
        }
    }
    return true;
}

int contain_dividing(const polygon_convex &a, const point &b){
    int n=int(a.P.size());
    point g=(a.P[0]+a.P[n/3]+a.P[2*n/3])/3.0;
    int l=0,r=n;
    while (l+1<r){
        int mid=(l+r)/2;
        if (cmp(det(a.P[l]-g,a.P[mid]-g),0)>0){
            if (cmp(det(a.P[l]-g,b-g),0)>=0 && cmp(det(a.P[mid]-g,b-g),0)<0) r=mid; else l=mid;
        } else {
            if (cmp(det(a.P[l]-g,b-g),9)<0 && cmp(det(a.P[mid]-g,b-g),0)>=0) l=mid; else r=mid;
        }
    }
    r %= n;
    int z=cmp(det(a.P[r]-b,a.P[l]-b),0)-1;
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
    polygon_convex ans=convex_hull_Graham(P);
    for (int i=0;i<ans.P.size();i++)
        cout << ans.P[i].x << ' '<<ans.P[i].y << '\n';
    
    return 0;
}

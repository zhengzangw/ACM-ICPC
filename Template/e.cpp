#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>
#define eps 1e-8
#define fabs(x) ((x) > 0 ? (x) : -(x))
#define zero(x) (fabs(x) < eps)
#define _sign(x) ((x) > eps ? 1 : ((x) < -eps) ? 2 : 0)
#define sqr(x) ((x) * (x))
const double pi = acos(-1);
struct point
{ //点、向量
    double x, y;
    int index;
    double ang;
    point() : x(0), y(0) {}
    point(double a, double b) : x(a), y(b){};
    void read() { scanf("%lf %lf", &x, &y); }
    friend point operator+(const point &a, const point &b)
    {
        return point(a.x + b.x, a.y + b.y);
    }
    friend point operator-(const point &a, const point &b)
    {
        return point(a.x - b.x, a.y - b.y);
    }
    bool operator==(const point &p) const
    {
        return zero(x - p.x) && zero(y - p.y);
    }
    friend point operator*(const point &a, const double &b)
    {
        return point(a.x * b, a.y * b);
    }
    friend point operator*(const double &a, const point &b)
    {
        return point(b.x * a, b.y * a);
    }
    friend point operator/(const point &a, const double &b)
    {
        return point(a.x / b, a.y / b);
    }
    bool operator<(const point &b) const
    {
        if (b.x == x)
            return y < b.y;
        return x < b.x;
    }
    void getangle()
    {
        ang = atan(y / x);
    }
};
typedef const point CP;
bool cmp(const point &p1, const point &p2)
{
    return p1.ang < p2.ang;
}

//外乘
double xmult(CP &p1, CP &p2, CP &p0)
{
    return (p1.x - p0.x) * (p2.y - p0.y) - (p2.x - p0.x) * (p1.y - p0.y);
}
double xmult(double x1, double y1, double x2, double y2, double x0, double y0)
{
    return (x1 - x0) * (y2 - y0) - (x2 - x0) * (y1 - y0);
}
double xmult(CP &v1, CP &v2)
{
    return v1.x * v2.y - v2.x * v1.y;
}

//内积
double dmult(CP &p1, CP &p2, CP &p0)
{
    return (p1.x - p0.x) * (p2.x - p0.x) + (p1.y - p0.y) * (p2.y - p0.y);
}
double dmult(double x1, double y1, double x2, double y2, double x0, double y0)
{
    return (x1 - x0) * (x2 - x0) + (y1 - y0) * (y2 - y0);
}
double dmult(CP &v1, CP &v2)
{
    return v1.x * v2.x + v1.y * v2.y;
}
struct line
{ //直线、线段
    double ang;
    point a, b;
    bool segment;
    line() {}
    line(const point &x, const point &y) : a(x), b(y), segment(1) {}
    line(const point &x, const point y, bool flag) : a(x), b(y), segment(flag) {}
    point line2point() { return point(b - a); }
    bool operator<(const line &y) const
    {
        if (zero(ang - y.ang))
            return (xmult(a, y.b, y.a) < 0);
        return ang < y.ang;
    }
    void getang()
    {
        ang = atan2(b.y - a.y, b.x - a.x);
    }
};
//距离
double dis(double x1, double y1, double x2, double y2)
{
    return sqrt(sqr(x1 - x2) + sqr(y1 - y2));
}
double dis(CP &p1, CP &p2)
{
    return sqrt(sqr(p1.x - p2.x) + sqr(p1.y - p2.y));
}
double dis2(double x1, double y1, double x2, double y2)
{
    return sqr(x1 - x2) + sqr(y1 - y2);
}
double dis2(CP &p1, CP &p2)
{
    return sqr(p1.x - p2.x) + sqr(p1.y - p2.y);
}
//点到直线距离
double disptoline(CP &p, CL &l)
{
    return fabs(xmult(p, l.a, l.b) / dis(l.a, l.b));
}
double disptoline(CP &p, CP &l1, CP &l2)
{
    return fabs(xmult(p, l1, l2) / dis(l1, l2));
}
typedef const line CL;
//圆
struct circle
{
    double r;
    point c;
    circle(){};
    circle(const point &p, double x)
    {
        c = p;
        r = x;
    }
};
typedef const circle CC;
//扇形面积
double p_circle_angle(CP &p1, CP &p2, CP &c, double r)
{
    double alpha = fabs(atan2(p1.y - c.y, p1.x - c.x) - atan2(p2.y - c.y, p2.x - c.x));
    return alpha * r;
}
//直线与圆
//判直线和圆是否有交点
int intersect_line_circle(CP &c, double r, CP &l1, CP &l2)
{
    return disptoline(c, l1, l2) < r + eps;
}
int intersect_line_circle(CC &c, CP &l1, CP &l2)
{
    return disptoline(c.c, l1, l2) < c.r + eps;
}
//圆与圆是否有交点
int intersect_circle_circle(CC &c1, CC &c2)
{
    return dis(c1.c, c2.c) < c1.r + c2.r + eps && dis(c1.c, c2.c) > fabs(c1.r - c2.r) - eps;
}
point intersection(CP &u1, CP &u2, CP &v1, CP &v2)
{
    point ret = u1;
    double t = ((u1.x - v1.x) * (v1.y - v2.y) - (u1.y - v1.y) * (v1.x - v2.x)) / ((u1.x - u2.x) * (v1.y - v2.y) - (u1.y - u2.y) * (v1.x - v2.x));
    ret.x += (u2.x - u1.x) * t;
    ret.y += (u2.y - u1.y) * t;
    return ret;
}
void intersection_line_circle(CP &c, double r, CP &l1, CP &l2, point &p1, point &p2)
{
    point p = c;
    double t;
    p.x += l1.y - l2.y;
    p.y += l2.x - l1.x;
    p = intersection(p, c, l1, l2);
    t = sqrt(r * r - dis2(p, c)) / dis(l1, l2);
    p1.x = p.x + (l2.x - l1.x) * t;
    p1.y = p.y + (l2.y - l1.y) * t;
    p2.x = p.x - (l2.x - l1.x) * t;
    p2.y = p.y - (l2.y - l1.y) * t;
}
void intersection_circle_circle(CC &c1, CC &c2, point &p1, point &p2)
{
    point u, v;
    double t;
    t = (1 + (sqr(c1.r) - sqr(c2.r)) / dis2(c1.c, c2.c)) / 2;
    u.x = c1.c.x + (c2.c.x - c1.c.x) * t;
    u.y = c1.c.y + (c2.c.y - c1.c.y) * t;
    v.x = u.x + c1.c.y - c2.c.y;
    v.y = u.y - c1.c.x + c2.c.x; //根轴line(u,v)
    intersection_line_circle(c1.c, c1.r, u, v, p1, p2);
}

circle c[1010];
int main(){
    int T,R,m,x,y,r;
    circle C;
    scanf("%d",&T);
    double ans;
    while (T--){
        scanf("%d%d",&m,&R);
        C.c = point(0, 0);
        C.r = R;
        ans = 2*pi*R;
        for (int i=1;i<=m;++i){
            scanf("%d%d%d",&x,&y,&r);
            c[i].c = point(x,y);
            c[i].r = r;
        }
        for (int i=1;i<=m;++i){
            if (!intersect_circle_circle(C,c[i])) continue;
            point p1,p2;
            intersection_circle_circle(C,c[i],p1,p2);
            ans += p_circle_angle(p1,p2,c[i].c,c[i].r);
            ans -= p_circle_angle(p1,p2,C.c,C.r);
        }
        printf("%.6f\n",ans)
    }
    return 0;
}
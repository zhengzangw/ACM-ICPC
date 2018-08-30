//Geometry.h
//仅供参考的计算几何模板
//Team NJU99
#ifndef Geometry_H
#define Geometry_H
#include <cmath>
#include <cstdio>
#include <iostream>
#include <algorithm>
#include <vector>
#include <complex>
#include <utility>
#define eps 1e-8
#define fabs(x) ((x) > 0 ? (x) : -(x))
#define zero(x) (fabs(x) < eps)
#define _sign(x) ((x) > eps ? 1 : ((x) < -eps) ? 2 : 0)
#define sqr(x) ((x) * (x))
const double pi = acos(-1);
//using namespace std;
#ifndef Geometry_2D_H
#define Geometry_2D_H

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

typedef const line CL;

//长度
double len(CP &v)
{
    return sqrt((sqr(v.x) + sqr(v.y)));
}
double len(CL &v)
{
    return sqrt(sqr(v.b.y - v.a.y) + sqr(v.b.x - v.a.x));
}

//三点共线
int dots_inline(CP &p1, CP &p2, CP &p3)
{
    return zero(xmult(p1, p2, p3));
}
int dots_inline(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return zero(xmult(x1, y1, x2, y2, x3, y3));
}
int dots_inline(CP &p, CL &l)
{
    return dots_inline(p, l.a, l.b);
}

//点在闭线段上
int dot_online_in(CP &p, CL &l)
{
    return zero(xmult(p, l.a, l.b)) && (l.a.x - p.x) * (l.b.x - p.x) < eps && (l.a.y - p.y) * (l.b.y) < eps;
}

int dot_online_in(CP &p, CP &l1, CP &l2)
{
    return zero(xmult(p, l1, l2)) && (l1.x - p.x) * (l2.x - p.x) < eps && (l1.y - p.y) * (l2.y - p.y) < eps;
}

int dot_online_in(double x, double y, double x1, double y1, double x2, double y2)
{
    return zero(xmult(x, y, x1, y1, x2, y2)) && (x1 - x) * (x2 - x) < eps && (y1 - y) * (y2 - y) < eps;
}

//点在开线段上
int dot_online_ex(CP &p, CL &l)
{
    return dot_online_in(p, l) && (!zero(p.x - l.a.x) || !zero(p.y - l.a.y)) && (!zero(p.x - l.b.x) || !zero(p.y - l.b.y));
}
int dot_online_ex(CP &p, CP &l1, CP &l2)
{
    return dot_online_in(p, l1, l2) && (!zero(p.x - l1.x) || !zero(p.y - l1.y)) && (!zero(p.x - l2.x) || !zero(p.y - l2.y));
}
int dot_online_ex(double x, double y, double x1, double y1, double x2, double y2)
{
    return dot_online_in(x, y, x1, y1, x2, y2) && (!zero(x - x1) || !zero(y - y1)) && (!zero(x - x2) || !zero(y - y2));
}

//点在同侧
int same_side(CP &p1, CP &p2, CL &l)
{
    return xmult(l.a, p1, l.b) * xmult(l.a, p2, l.b) > eps;
}
int same_side(CP &p1, CP &p2, CP &l1, CP &l2)
{
    return xmult(l1, p1, l2) * xmult(l1, p2, l2) > eps;
}
//点在异侧
int opposite_side(CP &p1, CP &p2, CL &l)
{
    return xmult(l.a, p1, l.b) * xmult(l.a, p2, l.b) < -eps;
}
int opposite_side(CP &p1, CP &p2, CP &l1, CP &l2)
{
    return xmult(l1, p1, l2) * xmult(l1, p2, l2) < -eps;
}
//平行
int parallel(CL &u, CL &v)
{
    return zero(xmult(u.a - u.b, v.a - v.b));
}
int parallel(CP &u1, CP &u2, CP &v1, CP &v2)
{
    return zero((u1.x - u2.x) * (v1.y - v2.y) - (v1.x - v2.x) * (u1.y - u2.y));
}
//垂直
int perpendicular(CL &u, CL &v)
{
    return zero((u.a.x - u.b.x) * (v.a.x - v.b.x) + (u.a.y - u.b.y) * (v.a.y - v.b.y));
}
int perpendicular(CP &u, CP &v)
{
    return zero((u.x * v.y + u.y * v.x));
}
//相交（包括端点和部分重合）
int intersect_in(CL &u, CL &v)
{
    if (!dots_inline(u.a, u.b, v.a) || !dots_inline(u.a, u.b, v.b))
        return !same_side(u.a, u.b, v) && !same_side(v.a, v.b, u);
    return dot_online_in(u.a, v) || dot_online_in(u.b, v) || dot_online_in(v.a, u) || dot_online_in(v.b, u);
}
//相交（不包括端点和部分重合）
int intersect_ex(CL &u, CL &v)
{
    return opposite_side(u.a, u.b, v) && opposite_side(v.a, v.b, u);
}
int intersect_ex(CP &u1, CP &u2, CP &v1, CP &v2)
{
    return opposite_side(u1, u2, v1, v2) && opposite_side(v1, v2, u1, u2);
}
//交点 (前提：相交)
//利用面积比计算相似比
point intersection(CL &u, CL &v)
{
    point ret = u.a;
    double t = t = ((u.a.x - v.a.x) * (v.a.y - v.b.y) - (u.a.y - v.a.y) * (v.a.x - v.b.x)) / ((u.a.x - u.b.x) * (v.a.y - v.b.y) - (u.a.y - u.b.y) * (v.a.x - v.b.x));
    ret.x += (u.b.x - u.a.x) * t;
    ret.y += (u.b.y - u.a.y) * t;
    return ret;
}
point intersection(CP &u1, CP &u2, CP &v1, CP &v2)
{
    point ret = u1;
    double t = ((u1.x - v1.x) * (v1.y - v2.y) - (u1.y - v1.y) * (v1.x - v2.x)) / ((u1.x - u2.x) * (v1.y - v2.y) - (u1.y - u2.y) * (v1.x - v2.x));
    ret.x += (u2.x - u1.x) * t;
    ret.y += (u2.y - u1.y) * t;
    return ret;
}

//中垂线
line pbline(CL &l)
{
    line ret;
    ret.a = (l.a + l.b) / 2;
    double a = l.b.x - l.a.x, b = l.b.y - l.a.y;
    double c = (l.a.y - l.b.y) * ret.a.y + (l.a.x - l.b.x) * ret.a.x;
    if (!zero(a))
    {
        ret.b.y = 0;
        ret.b.x = -c / a;
        if (zero(dis(ret.a, ret.b)))
        {
            ret.b.y = 1e10;
            ret.b.x = -(c - b * ret.b.y) / a;
        }
    }
    else
    {
        ret.b.x = 0.0;
        ret.b.y = -c / b;
        if (zero(dis(ret.a, ret.b)))
        {
            ret.b.x = 1e10;
            ret.b.y = -(c - a * ret.b.x) / b;
        }
    }
    return ret;
}
//角平分线
line jpfline(CP &a, CP &b, CP &c)
{
    point ua, ub;
    double m, n;
    ua = a;
    m = atan2(b.y - a.y, b.x - a.x);
    n = atan2(c.y - a.y, c.x - a.x);
    ub.x = ua.x + cos((m + n) / 2);
    ub.y = ua.y + sin((m + n) / 2);
    return line(ua, ub);
}

//垂足
point ptoline(CP &p, CL &l)
{
    point t = p;
    t.x += l.a.y - l.b.y;
    t.y += l.b.x - l.a.x;
    return intersection(p, t, l.a, l.b);
}
point ptoline(CP &p, CP &l1, CP &l2)
{
    point t = p;
    t.x += l1.y - l2.y;
    t.y += l1.x - l2.x;
    return intersection(p, t, l1, l2);
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
//点到直线最近点
point ptoseg(CP &p, CL &l)
{
    point t = p;
    t.x += l.a.y - l.b.y;
    t.y += l.b.x - l.a.x;
    if (xmult(l.a, t, p) * xmult(l.b, t, p) > eps)
        return dis(p, l.a) < dis(p, l.b) ? l.a : l.b;
    return intersection(p, t, l.a, l.b);
}
point ptoseg(CP &p, CP &l1, CP &l2)
{
    point t = p;
    t.x += l1.y - l2.y;
    t.y += l1.x - l2.x;
    if (xmult(l1, t, p) * xmult(l2, t, p) > eps)
        return dis(p, l1) < dis(p, l2) ? l1 : l2;
    return intersection(p, t, l1, l2);
}
//点到线段距离
double disptoseg(CP &p, CL &l)
{
    point t = p;
    t.x += l.a.y - l.b.y;
    t.y += l.b.x - l.a.x;
    if (xmult(l.a, t, p) * xmult(l.b, t, p) > eps)
        return dis(p, l.a) < dis(p, l.b) ? dis(p, l.a) : dis(p, l.b);
    return fabs(xmult(p, l.a, l.b) / dis(l.a, l.b));
}
double disptoseg(CP &p, CP &l1, CP &l2)
{
    point t = p;
    t.x += l1.y - l2.y;
    t.y += l2.x - l1.x;
    if (xmult(l1, t, p) * xmult(l2, t, p) > eps)
        return dis(p, l1) < dis(p, l2) ? dis(p, l1) : dis(p, l2);
    return fabs(xmult(p, l1, l2) / dis(l1, l2));
}
//线段到线段距离
double dissegtoseg(CL &l1, CL &l2)
{
    if (intersect_in(l1, l2))
        return 0;
    return fmin(fmin(disptoseg(l1.a, l2), disptoseg(l1.b, l2)), fmin(disptoseg(l2.a, l1), disptoseg(l2.b, l1)));
}
double dissegtoseg(CP &l1a, CP &l1b, CP &l2a, CP &l2b)
{
    if (intersect_in(line(l1a, l1b), line(l2a, l2b)))
        return 0;
    return fmin(fmin(disptoseg(l1a, l2a, l2b), disptoseg(l1b, l2a, l2b)), fmin(disptoseg(l2a, l1a, l1b), disptoseg(l2b, l1a, l1b)));
}
//旋转
//矢量V 以P 为顶点逆时针旋转angle 并放大scale 倍
point rotate(point v, point p, double angle, double scale)
{
    point ret = p;
    v.x -= p.x;
    v.y -= p.y;
    p.x = scale * cos(angle);
    p.y = scale * sin(angle);
    ret.x += v.x * p.x - v.y * p.y;
    ret.y += v.x * p.y + v.y * p.x;
    return ret;
}

//坐标变换 在新坐标系O(I,e1,e2) 中的坐标
point rotate(CP &p, CP &I, CP &e1, CP &e2)
{
    point p2;
    p2.x = I.x + e1.x * p.x + e1.y * p.y;
    p2.y = I.y + e2.x * p.x + e2.y * p.y;
    return p2;
}

point rotate(CP &p, double angle)
{
    point e1, e2, I;
    e1.x = cos(angle);
    e1.y = -sin(angle);
    e2.x = -e1.y;
    e2.y = e1.x;
    I.x = 0;
    I.y = 0;
    return rotate(p, I, e1, e2);
}

//旋转角
double angle(CP &v1, CP &v2)
{
    double cosa = dmult(v1, v2) / len(v1) / len(v2);
    cosa = acos(cosa);
    if (cosa*2>pi)
        cosa = 2*pi-cosa;
    return cosa;
}
double angle(CP &a, CP &b, CP &c)
{
    return angle(a - b, b - c);
}

//外心
point circumcenter(CP &a, CP &b, CP &c)
{
    point ua, ub, va, vb;
    ua.x = (a.x + b.x) / 2;
    ua.y = (a.y + b.y) / 2;
    ub.x = ua.x - a.y + b.y;
    ub.y = ua.y + a.x - b.x;
    va.x = (a.y + c.y) / 2;
    vb.x = va.x - a.y + c.y;
    vb.y = va.y + a.x - c.x;
    return intersection(ua, ub, va, vb);
}
//内心
point incenter(CP &a, CP &b, CP &c)
{
    point ua, ub, va, vb;
    double m, n;
    ua = a;
    m = atan2(b.y - a.y, b.x - a.x);
    n = atan2(c.y - a.y, c.x - a.x);
    ub.x = ua.x + cos((m + n) / 2);
    ub.y = ua.y + sin((m + n) / 2);
    va = b;
    m = atan2(a.y - b.y, a.x - b.x);
    n = atan2(c.y - b.y, c.x - b.x);
    vb.x = va.x + cos((m + n) / 2);
    vb.y = va.y + sin((m + n) / 2);
    return intersection(ua, ub, va, vb);
}
//垂心
point perpencenter(CP &a, CP &b, CP &c)
{
    point ua, ub, va, vb;
    ua = c;
    ub.x = ua.x - a.y + b.y;
    ub.y = ua.y + a.x - b.x;
    va = b;
    vb.x = va.x - a.y + c.y;
    vb.y = va.y + a.x - c.x;
    return intersection(ua, ub, va, vb);
}
//重心（三边距离平方和最小，距离积最大）
point barycenter(CP &a, CP &b, CP &c)
{
    point ua, ub, va, vb;
    ua.x = (a.x + b.x) / 2;
    ua.y = (a.y + b.y) / 2;
    ub = c;
    va.x = (a.x + c.x) / 2;
    va.y = (a.y + c.y) / 2;
    vb = b;
    return intersection(ua, ub, va, vb);
}
//费马点/等角中心
//到三角形三顶点距离之和最小的点
/*若有一个角大于120度。那么这个角所在的点就是费马点。
 若不存在。那么对于三角形ABC，任取两条边（如果AB、AC），向外做等边三角形得到C'和 A'。那么AA'和CC'的交点就是费马点。 */
point fermatpoint(CP &a, CP &b, CP &c)
{
    //模拟退火
    point u, v;
    double step = fabs(a.x) + fabs(a.y) + fabs(b.x) + fabs(b.y) + fabs(c.x) + fabs(c.y);
    int i, j, k;
    u = (a + b + c) / 3;
    while (step > 1e-10)
        for (k = 0; k < 10; step /= 2, k++)
            for (i = -1; i <= 1; i++)
                for (j = -1; j <= 1; j++)
                {
                    v.x = u.x + step * i;
                    v.y = u.y + step * j;
                    if (dis(u, a) + dis(u, b) + dis(u, c) > dis(v, a) + dis(v, b) + dis(v, c))
                        u = v;
                }
    return u;
}

//判定凸多边形
//顶点按顺时针或逆时针给 出, 允许相邻边共线
int is_convex(int n, point *p)
{
    int i, s[3] = {1, 1, 1};
    for (i = 0; i < n && s[1] | s[2]; i++)
        s[_sign(xmult(p[(i + 1) % n], p[(i + 2) % n], p[i]))] = 0;
    return s[1] | s[2];
}
//不允许相邻边共线
int is_convex2(int n, point *p)
{
    int i, s[3] = {1, 1, 1};
    for (i = 0; i < n && s[0] && s[1] | s[2]; i++)
        s[_sign(xmult(p[(i + 1) % n], p[(i + 2) % n], p[i]))] = 0;
    return s[0] && s[1] | s[2];
}
//点与多边形关系
//判点在凸多边形内或多边形边上, 顶点按顺时 针或逆时针给出
int inside_convex(CP &q, int n, point *p)
{
    int i, s[3] = {1, 1, 1};
    for (i = 0; i < n && s[1] | s[2]; i++)
        s[_sign(xmult(p[(i + 1) % n], q, p[i]))] = 0;
    return s[1] | s[2];
}
//判点在凸多边形内, 顶点按顺时针或逆时针给出, 在多边形边上返回0
int inside_convex2(CP &q, int n, point *p)
{
    int i, s[3] = {1, 1, 1};
    for (i = 0; i < n && s[0] && s[1] | s[2]; i++)
        s[_sign(xmult(p[(i + 1) % n], q, p[i]))] = 0;
    return s[0] && s[1] | s[2];
}
//判点在任意多边形内顶点按顺时针或逆时针给出,
//on_edge表示点在多边形边上时的返回值,offset为多边形坐标上限
int offset;
int inside_polygon(CP &q, int n, point *p, int on_edge = 1)
{
    //射线法
    point q2;
    int i = 0, count = 0;
    while (i < n)
        for (count = i = 0, q2.x = rand() + offset, q2.y = rand() + offset; i < n; ++i)
            if (zero(xmult(q, p[i], p[(i + 1) % n])) && (p[i].x - q.x) * (p[(i + 1) % n].x - q.x) < eps && (p[i].y - q.y) * (p[(i + 1) % n].y - q.y) < eps)
                return on_edge;
            else if (zero(xmult(q, q2, p[i])))
                break;
            else if (xmult(q, p[i], q2) * xmult(q, p[(i + 1) % n], q2) < -eps && xmult(p[i], q, p[(i + 1) % n]) * xmult(p[i], q2, p[(i + 1) % n]) < -eps)
                count++;
    return count & 1;
} //还有 面积法、角度法

//线段与多边形关系
//判线段在任意多边形内, 顶点按顺时针或逆时 427e 针给出, 与边界相交返回1
int MAXN;
int inside_polygonn(CP &l1, CP &l2, int n, point *p)
{
    point t[MAXN], tt;
    int i, j, k = 0;
    if (!inside_polygon(l1, n, p) || !inside_polygon(l2, n, p))
        return 0;
    for (i = 0; i < n; i++)
        if (opposite_side(l1, l2, p[i], p[(i + 1) % n]) && opposite_side(p[i], p[(i + 1) % n], l1, l2))
            return 0;
        else if (dot_online_in(l1, p[i], p[(i + 1) % n]))
            t[k++] = l1;
        else if (dot_online_in(l2, p[i], p[(i + 1) % n]))
            t[k++] = l2;
        else if (dot_online_in(p[i], l1, l2))
            t[k++] = p[i];
    for (i = 0; i < k; i++)
        for (j = i + 1; j < k; j++)
        {
            tt = (t[i] + t[j]) / 2;
            if (!inside_polygon(tt, n, p))
                return 0;
        }
    return 1;
}

//多边形重心
point barycenter(int n, point *p)
{
    point ret, t;
    double t1 = 0, t2;
    int i;
    ret.x = ret.y = 0;
    for (i = 1; i < n - 1; i++)
        if (fabs(t2 = xmult(p[0], p[i], p[i + 1])) > eps)
        {
            t = barycenter(p[0], p[i], p[i + 1]);
            ret.x += t.x * t2;
            ret.y += t.y * t2;
            t1 += 2;
        }
    if (fabs(t1) > eps)
        ret.x /= t1;
    ret.y /= t1;
    return ret;
}
//切割
//将多边形沿l1,l2确定的直线切割在side侧切割,保证l1,l2,side不共线
void polygon_cut(int &n, point *p, CP &l1, CP &l2, CP &side)
{
    point pp[100];
    int m = 0, i;
    for (i = 0; i < n; i++)
    {
        if (same_side(p[i], side, l1, l2))
            pp[m++] = p[i];
        if (!same_side(p[i], p[(i + 1) % n], l1, l2) && !(zero(xmult(p[i], l1, l2)) && zero(xmult(p[(i + 1) % n], l1, l2))))
            pp[m++] = intersection(p[i], p[(i + 1) % n], l1, l2);
    }
    for (n = i = 0; i < m; i++)
        if (!i || !zero(pp[i].x - pp[i - 1].x) || !zero(pp[i].y - pp[i - 1].y))
            p[n++] = pp[i];
    if (zero(p[n - 1].x - p[0].x) && zero(p[n - 1].y - p[0].y))
        n--;
    if (n < 3)
        n = 0;
}
//最长线段

//三角形面积
double area_triangle(CP &p1, CP &p2, CP &p3)
{
    return xmult(p1, p2, p3) / 2;
}
double area(CP &p1, CP &p2, CP &p3)
{
    return fabs(xmult(p1, p2, p3)) / 2;
}
double area(double x1, double y1, double x2, double y2, double x3, double y3)
{
    return fabs(xmult(x1, y1, x2, y2, x3, y3)) / 2;
}
double area(double a, double b, double c)
{
    //海伦公式
    double s = (a + b + c) / 2;
    return sqrt(s * (s - a) * (s - b) * (s - c));
}
//多边形面积
//顶点按顺时针或逆时针给出,逆为正面积
double area_polygon(int n, point *p)
{
    double s1 = 0, s2 = 0;
    int i;
    for (i = 0; i < n; i++)
    {
        s1 += p[(i + 1) % n].y * p[i].x;
        s2 += p[(i + 1) % n].y * p[(i + 2) % n].x;
    }
    return (s1 - s2) / 2;
}

//求平行于v的所有射线中，穿过的凸包中最左边的点的坐标;凸包点按顺时针给出
point vector_throw_convex(int n, point *convex, CP &v)
{
    int s = 0;
    double as = angle(v, convex[s], convex[(s + 1) % n]);
    int t = n - 1;
    double at = angle(v, convex[t], convex[(t + 1) % n]);
    while (s < t)
    {
        if (as >= at)
        {
            s = t;
            break;
        }
        int mid = (s + t + 1) / 2;
        double amid = angle(v, convex[mid], convex[(mid + 1) % n]);
        if (amid <= as)
        {
            s = mid;
            as = amid;
        }
        else
        {
            t = mid - 1;
            at = angle(v, convex[t], convex[(t + 1) % n]);
        }
    }
    return convex[(s + 1) % n];
}

//求直线l1是否穿过凸包, 凸包按顺时针给出，返回是否穿过
// p储存凸包内与l1共线的某点
bool line_throw_convex(int n, point *convex, CL &l, point &p)
{
    point p1 = vector_throw_convex(n, convex, l.a - l.b);
    point p2 = vector_throw_convex(n, convex, l.b - l.a);
    line l2(p1, p2);
    p = intersection(l, l2);
    if (dot_online_in(p, l2))
        return true;
    return false;
}
//求射线是否穿过凸包, 凸包按顺时针给出，返回是否穿过
bool ray_throw_convex(int n, point *convex, CL &l, point &p)
{
    if (line_throw_convex(n, convex, l, p))
        if (dmult(p, l.b, l.a) >= -eps)
            return true;
    return false;
}

//凸包直径
//顺时针输入凸包，没有共线的点
double convex_diameter(int n, point *convex)
{
    //rotate_calipers
    int q = 1;
    double ans = 0;
    for (int p = 0; p < n; ++p)
    {
        while (xmult(convex[(p + 1) % n], convex[(q + 1) % n], convex[p]) < xmult(convex[(p + 1) % n], convex[q], convex[p]))
            q = (q + 1) % n;
        ans = fmax(ans, fmax(dis(convex[p], convex[q]), dis(convex[(p + 1) % n], convex[(q + 1) % n])));
    }
    return ans;
}
//凸包最小截面
//顺时针输入凸包
double convex_min_section(int n, point *convex)
{
    double l1 = 1000000000;
    for (int i = 0; i < n; ++i)
    {
        point a = convex[i] - convex[(i + 1) % n];
        point b = vector_throw_convex(n, convex, a);
        l1 = fmin(l1, disptoline(b, convex[i], convex[i + 1]));
    }
    return l1;
}
//凸包距离,逆时针输入
double convex_min_dis(int n, point *a, int m, point *b)
{
    //旋转卡壳
    int p1 = 0, p2 = 0;
    double ans = 1 << 30;
    for (int i = 0; i < n; ++i)
        if (a[i].y < a[p1].y)
            p1 = i;
    for (int i = 0; i < m; ++i)
        if (b[i].y > b[p1].y)
            p2 = i;
    for (int i = 0; i < n; ++i)
    {
        double t = xmult(b[(p2 + 1) % m], a[p1], a[(p1 + 1) % n]);
        t -= xmult(b[p2], a[p1], a[(p1 + 1) % n]);
        if (_sign(t) == 1)
        {
            ans = fmin(ans, disptoseg(a[p1], b[p2], b[(p2 + 1) % m]));
            p2 = (p2 + 1) % m;
            --i;
        }
        else if (_sign(t) == 2)
        {
            ans = fmin(ans, disptoseg(b[p2], a[p1], a[(p1 + 1) % n]));
            p1 = (p1 + 1) % n;
        }
        else
        {
            ans = fmin(ans, dissegtoseg(a[p1], a[(p1 + 1) % n], b[p2], b[(p2 + 1) % m]));
            p1 = (p1 + 1) % n;
            p2 = (p2 + 2) % m;
        }
    }
    return ans;
}
//多边形中最长的线段的长度，线段储存在l中
double inside_polygon_max(int n, point *p, line &l)
{
    double len = 0;
    for (int i = 0; i < n; ++i)
        for (int j = i + 1; j < n; ++j)
        {
            std::vector<point> points;
            points.clear();
            points.push_back(p[i]);
            points.push_back(p[j]);
            for (int a = 0; a < n; ++a)
                for (int b = a + 1; b < n; ++b)
                {
                    if (a == i)
                        continue;
                    if (parallel(p[i], p[j], p[a], p[b]))
                        continue;
                    point p1 = intersection(p[i], p[j], p[a], p[b]);
                    if (dmult(p[a], p[b], p1) <= 0)
                        points.push_back(p1);
                }
            sort(points.begin(), points.end());

            int s = 0;
            for (int k = 0; k < points.size() - 1; ++k)
            {
                if (zero(dis(points[k], points[k + 1])))
                    continue;
                point p1;
                p1 = (points[k] + points[k + 1]) / 2;
                if (inside_polygon(p1, n, p))
                    continue;
                double d = dis(points[s], points[k + 1]);
                if (len < d)
                {
                    len = d;
                    l.a = points[s];
                    l.b = points[k + 1];
                }
                s = k + 1;
            }
            double d = dis(points[s], points[points.size() - 1]);
            if (len < d)
            {
                len = d;
                l.a = points[s];
                l.b = points[points.size() - 1];
            }
        }

    return len;
}

//凸包生成graham O(nlogn)
point p1, p2;
int graham_cp(const void *a, const void *b)
{
    double ret = xmult(*((point *)a), *((point *)b), p1);
    return zero(ret) ? (xmult(*((point *)a), *((point *)b), p2) > 0 ? 1 : -1) : (ret > 0 ? 1 : -1);
}
void _graham(int n, point *p, int &s, point *ch)
{
    int i, k = 0;
    for (p1 = p2 = p[0], i = 1; i < n; p2.x += p[i].x, p2.y += p[i].y, i++)
        if (p1.y - p[i].y > eps || (zero(p1.y - p[i].y) && p1.x > p[i].x))
            p1 = p[k = i];
    p2.x /= n;
    p2.y /= n;
    p[k] = p[0];
    p[0] = p1;
    qsort(p + 1, n - 1, sizeof(point), graham_cp);
    for (ch[0] = p[0], ch[1] = p[1], ch[2] = p[2], s = i = 3; i < n; ch[s++] = p[i++])
        for (; s > 2 && xmult(ch[s - 2], p[i], ch[s - 1]) < -eps; s--)
            ;
}

//生成整理好的凸包
//返回凸包大小,凸包的点在convex中。
//参数maxsize为1则包含共线点,为0则不包含共线点,缺省为1;
//参数clockwise为1则顺时针构造,为0则逆时针构造,缺省为1
//在输入仅有若干共线点时算法不稳定,可能有此类情况请另行处理!;不能去掉点集中重合的点
int graham(int n, point *p, point *convex, int maxsize = 1, int dir = 1)
{
    point *temp = new point[n];
    int s, i;
    _graham(n, p, s, temp);
    for (convex[0] = temp[0], n = 1, i = (dir ? 1 : (s - 1)); dir ? (i < s) : i; i += (dir ? 1 : -1))
        if (maxsize || !zero(xmult(temp[i - 1], temp[i], temp[(i + 1) % s])))
            convex[n++] = temp[i];
    delete[] temp;
    return n;
}

//点与半平面
//判断点知否在半平面内,平面位于向量左侧
bool phplanerout(CL &l, CP &p)
{
    return xmult(p, l.b, l.a) > eps;
}
//半平面交，平面位于向量左侧 O(NlogN)
//另有分治发O(NlogN)
//多边形的核：每条边延长后，作半平面交
int halfpanelcross(int n, line *lines, point *p)
{
    int i;
    for (i = 0; i < n; ++i)
        lines[i].getang();
    std::sort(lines, lines + n);

    int m = 1;
    for (i = 1; i < n; ++i)
        if (!zero(lines[i].ang - lines[i - 1].ang))
            lines[m++] = lines[i];
    n = m;

    int bot = 0, top = 1;
    for (i = 2; i < n; ++i)
    {
        if ((parallel(lines[top], lines[top - 1]) || parallel(lines[bot], lines[bot + 1])))
            return 0;
        while ((bot < top) && (xmult(intersection(lines[top], lines[top - 1]), lines[i].b, lines[i].a) > eps))
            --top;
        while ((bot < top) && (xmult(intersection(lines[bot], lines[bot + 1]), lines[i].b, lines[i].a) > eps))
            ++bot;
        ++top;
        lines[top] = lines[i];
    }
    while ((bot < top) && (xmult(intersection(lines[top], lines[top - 1]), lines[bot].b, lines[bot].a) > eps))
        --top;
    while ((bot < top) && (xmult(intersection(lines[bot], lines[bot + 1]), lines[top].b, lines[top].a) > eps))
        ++bot;

    if (top <= bot + 1)
        return 0;
    n = 0;
    for (int i = bot; i < top; ++i)
        p[n++] = intersection(lines[i], lines[i + 1]);
    if (bot < top + 1)
        p[n++] = intersection(lines[bot], lines[top]);
    return n;
}

//折射
//u关于边v的折射角，不考虑全反射，折射率为r
line refraction(CL &u, CL &v, CP &p, double r)
{
    line v2;
    v2.a = p;
    v2.b.x = v.b.y - v.a.y + p.x;
    v2.b.y = v.a.x - v.b.x + p.x;
    if (dmult(v2.b - v2.a, u.b - u.a) < 0)
        std::swap(v2.a, v2.b);

    double alpha = xmult(v2.b - v2.a, u.b - u.a) / len(v2.b - v2.a) / len(u.b - u.a);
    alpha = asin(alpha / r) + atan2(v2.b.y - v2.a.y, v2.b.x - v2.a.x);
    v2.a = p;
    v2.b.x = 10 * cos(alpha) + v2.a.x;
    v2.b.y = 10 * sin(alpha) + v2.a.y;
    return v2;
}

//v关于凸包折射两次的情况，出射角保存在v中，折射率为r;如果不相交,返回false。不考虑镜面反射和入射到凸包角上的情况
bool refraction(int n, point p[], line &v, double r)
{
    int index = -1;
    line l[n];
    point p1;
    for (int i = 0; i < n; ++i)
    {
        l[i].a = p[i];
        l[i].b = p[(i + 1) % n];
    }
    for (int i = 0; i < n; ++i)
    {
        if (parallel(v, l[i]))
            continue;
        point p2 = intersection(v, l[i]);
        if (dmult(l[i].a, l[i].b, p2) >= 0 || dmult(p2, v.b, v.a) <= 0)
            continue;
        if ((index == -1) || (dis2(p2, v.a) < dis2(p1, v.a)))
        {
            index = i;
            p1 = p2;
        }
    }
    if (index == -1)
        return 0;
    std::swap(l[0], l[index]);
    v = refraction(v, l[0], p1, r);
    index = -1;
    for (int i = 0; i < n; ++i)
    {
        if (parallel(v, l[i]))
            continue;
        point p2 = intersection(v, l[i]);
        if (dmult(l[i].a, l[i].b, p2) >= 0 || dmult(p2, v.b, v.a) <= 0)
            continue;
        if ((index == -1) || (dis2(p2, v.a) < dis2(p1, v.a)))
        {
            index = i;
            p1 = p2;
        }
    }
    std::swap(l[1], l[index]);
    v = refraction(v, l[1], p1, 1 / r);
    return true;
}

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
//判线段和圆是否有交点
int intersect_seg_circle(CP &c, double r, CP &l1, CP &l2)
{
    double t1 = dis(c, l1) - r, t2 = dis(c, l2) - r;
    point t = c;
    if (t1 < eps || t2 < eps)
        return t1 > -eps || t2 > -eps;
    t.x += l1.y - l2.y;
    t.y += l2.x - l1.x;
    return xmult(l1, c, t) * xmult(l2, c, t) < eps && disptoline(c, l1, l2) - r < eps;
}
int intersect_seg_circle(CC &c, CP &l1, CP &l2)
{
    double t1 = dis(c.c, l1) - c.r, t2 = dis(c.c, l2) - c.r;
    point t = c.c;
    if (t1 < eps || t2 < eps)
        return t1 > -eps || t2 > -eps;
    t.x += l1.y - l2.y;
    t.y += l2.x - l1.x;
    return xmult(l1, c.c, t) * xmult(l2, c.c, t) < eps && disptoline(c.c, l1, l2) - c.r < eps;
}
//圆与圆是否有交点
int intersect_circle_circle(CP &c1, double r1, CP &c2, double r2)
{
    return dis(c1, c2) < r1 + r2 + eps && dis(c1, c2) > fabs(r1 - r2) - eps;
}
int intersect_circle_circle(CC &c1, CC &c2)
{
    return dis(c1.c, c2.c) < c1.r + c2.r + eps && dis(c1.c, c2.c) > fabs(c1.r - c2.r) - eps;
}
//计算圆上到点p最近点, 如p与圆心重合, 返回p本身
point dot_to_circle(CP &c, double r, CP &p)
{
    point u, v;
    if (dis(p, c) < eps)
        return p;
    u.x = c.x + r * fabs(c.x - p.x) / dis(c, p);
    u.y = c.y + r * fabs(c.y - p.y) / dis(c, p) * ((c.x - p.x) * (c.y - p.y) < 0 ? -1 : 1);
    v.x = c.x - r * fabs(c.x - p.x) / dis(c, p);
    v.y = c.y - r * fabs(c.y - p.y) / dis(c, p) * ((c.x - p.x) * (c.y - p.y) < 0 ? -1 : 1);
    return dis(u, p) < dis(v, p) ? u : v;
}
point dot_to_circle(CC &c, CP &p)
{
    point u, v;
    if (dis(p, c.c) < eps)
        return p;
    u.x = c.c.x + c.r * fabs(c.c.x - p.x) / dis(c.c, p);
    u.y = c.c.y + c.r * fabs(c.c.y - p.y) / dis(c.c, p) * ((c.c.x - p.x) * (c.c.y - p.y) < 0 ? -1 : 1);
    v.x = c.c.x - c.r * fabs(c.c.x - p.x) / dis(c.c, p);
    v.y = c.c.y - c.r * fabs(c.c.y - p.y) / dis(c.c, p) * ((c.c.x - p.x) * (c.c.y - p.y) < 0 ? -1 : 1);
    return dis(u, p) < dis(v, p) ? u : v;
}
//直线与圆交点
//保证直线与圆有交点;线段与圆的交点可用这个函数后判点是否在线段上
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
void intersection_line_circle(CC &c, CP &l1, CP &l2, point &p1, point &p2)
{
    point p = c.c;
    double t;
    p.x += l1.y - l2.y;
    p.y += l2.x - l1.x;
    p = intersection(p, c.c, l1, l2);
    t = sqrt(c.r * c.r - dis2(p, c.c)) / dis(l1, l2);
    p1.x = p.x + (l2.x - l1.x) * t;
    p1.y = p.y + (l2.y - l1.y) * t;
    p2.x = p.x - (l2.x - l1.x) * t;
    p2.y = p.y - (l2.y - l1.y) * t;
}

//点对圆的幂
double mi(CC &c, CP &p)
{
    return dis2(c.c, p) - c.r * c.r;
}
double mi(CP &c, double r, CP &p)
{
    return dis2(c, p) - r * r;
}
//圆与圆交点
//保证圆心不重合
void intersection_circle_circle(CP &c1, double r1, CP &c2, double r2, point &p1, point &p2)
{
    point u, v;
    double t;
    t = (1 + (sqr(r1) - sqr(r2)) / dis2(c1, c2)) / 2;
    u.x = c1.x + (c2.x - c1.x) * t;
    u.y = c1.y + (c2.y - c1.y) * t;
    v.x = u.x + c1.y - c2.y;
    v.y = u.y - c1.x + c2.x; //根轴line(u,v)
    intersection_line_circle(c1, r1, u, v, p1, p2);
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
//圆在多边形内
//顶点按顺时针或逆时针给出, offset为多边形坐标上限
bool inside_circle_polygon(CP &c, double r, int n, point *polygon)
{
    if (!inside_polygon(c, n, polygon, 1))
        return false;
    for (int i = 0; i < n; ++i)
        if (disptoline(c, polygon[i], polygon[(i + 1) % n]) < r)
            return false;
    return true;
}
bool inside_circle_polygon(CC &c, int n, point *polygon)
{
    if (!inside_polygon(c.c, n, polygon, 1))
        return false;
    for (int i = 0; i < n; ++i)
        if (disptoline(c.c, polygon[i], polygon[(i + 1) % n]) < c.r)
            return false;
    return true;
}
//多边形在圆内或上
bool inside_polygon_circle(CP &c, double r, int n, point *polygon)
{
    for (int i = 0; i < n; ++i)
        if (dis2(c, polygon[i]) >= r * r)
            return false;
    return true;
}
bool inside_polygon_circle(CC &c, int n, point *polygon)
{
    for (int i = 0; i < n; ++i)
        if (dis2(c.c, polygon[i]) >= c.r * c.r)
            return false;
    return true;
}
//求圆外一点到圆切线,返回两个切点
void tangent_point_circle(CP &c, double r, CP &p, point &a, point &b)
{
    double d = dis(c, p);
    double angp = acos(r / d);
    double ango = atan2(p.y - c.y, p.x - c.x);
    a.x = c.x + r * cos(ango + angp);
    a.y = c.y + r * sin(ango + angp);
    b.x = c.x + r * cos(ango - angp);
    b.y = c.y + r * sin(ango - angp);
}
void tangent_point_circle(CC &c, CP &p, point &a, point &b)
{
    double d = dis(c.c, p);
    double angp = acos(c.r / d);
    double ango = atan2(p.y - c.c.y, p.x - c.c.x);
    a.x = c.c.x + c.r * cos(ango + angp);
    a.y = c.c.y + c.r * sin(ango + angp);
    b.x = c.c.x + c.r * cos(ango - angp);
    b.y = c.c.y + c.r * sin(ango - angp);
}
//求两个圆的内公切线 (保证相离)
//同法苛求两个圆的外公切线，但有精度问题
void incut_circle_circle(CP &c1, double r1, CP &c2, double r2, line &l1, line &l2)
{
    double d = sqrt(dis2(c1, c2) - sqr(r1 + r2)); //内公切线长
    point p1, p2;
    intersection_circle_circle(c1, r1 + r2, c2, d, p1, p2);
    l1.a = (p1 * r1 + c1 * r2) / (r1 + r2);
    l1.b = l1.a + (c2 - p1);
    l2.a = (p2 * r1 + c1 * r2) / (r1 + r2);
    l2.b = l2.a + (c2 - p2);
}

//三角形外接圆
void out_circleoftri(CP &a, CP &b, CP &c, circle &tmp)
{
    tmp.c = circumcenter(a, b, c);
    tmp.r = dis(a, tmp.c);
}
//三角形内接圆
void in_circleoftri(CP &a, CP &b, CP &c, circle &tmp)
{
    tmp.c = incenter(a, b, c);
    tmp.r = disptoline(tmp.c, a, b);
}

//包含n个点的最小圆 (n<=3)
void min_circle_reduce(int n, point *p, circle &tmp)
{
    if (n == 0)
        tmp.r = -2;
    else if (n == 1)
    {
        tmp.c = p[0];
        tmp.r = 0;
    }
    else if (n == 2)
    {
        tmp.r = dis(p[0], p[1]) / 2;
        tmp.c = (p[0] + p[1]) / 2;
    }
    else if (n == 3)
        out_circleoftri(p[0], p[1], p[2], tmp);
}

void min_circle(int n, point *p, int m, point *down, circle &c)
{
    min_circle_reduce(m, down, c);
    if (m == 3)
        return;
    for (int i = 0; i < n; ++i)
    {
        if (dis(p[i], c.c) > c.r)
        {
            down[m] = p[i];
            min_circle(i, p, m + 1, down, c);
            point tmp = p[i];
            for (int j = i; j >= i; --j)
                p[j] = p[j - 1];
            p[0] = tmp;
        }
    }
}
//包含n个点的最小圆 O(n)
//p随机
void min_circle(int n, point *p, circle &c)
{
    point down[3];
    min_circle(n, p, 0, down, c);
}

//扇形面积
double area_circle_angle(CP &p1, CP &p2, CP &c, double r)
{
    double alpha = fabs(atan2(p1.y - c.y, p1.x - c.x) - atan2(p2.y - c.y, p2.x - c.x));
    if (alpha > pi)
        alpha = 2 * pi - alpha;
    return alpha / 2 * r * r;
}
//圆与三角形相交面积 （三角形p1p2c）
double area_triangle_circle(CP &c, double r, CP &p1, CP &p2)
{
    double x = xmult(p2, c, p1);
    int flag = _sign(x);
    if (flag == 0)
        return 0;
    double r2 = sqr(r);
    double s = 0, l1 = dis2(p1, c), l2 = dis2(p2, c);
    if ((l1 <= r2) && (l2 <= r2))
        return area_triangle(p2, c, p1) * flag;
    if ((l1 > r2) && (l2 > r2))
    {
        point p3, p4;
        s = area_circle_angle(p2, p1, c, r);
        if (disptoseg(c, p1, p2) < r)
        {
            intersection_line_circle(c, r, p1, p2, p3, p4);
            if (dis2(p3, p1) > dis2(p4, p1))
                std::swap(p3, p4);
            s -= area_circle_angle(p3, p4, c, r) - area_triangle(p3, c, p4);
        }
        return s * flag;
    }
    if (l1 < l2)
    {
        point p3, p4;
        intersection_line_circle(c, r, p1, p2, p3, p4);
        if (dmult(p3, p2, p1) <= 0)
            p3 = p4;
        s = area_triangle(p3, p2, c) + area_circle_angle(p3, p2, c, r);
        return s * flag;
    }
    else
    {
        point p3, p4;
        intersection_line_circle(c, r, p1, p2, p3, p4);
        if (dmult(p3, p2, p1) <= 0)
            p3 = p4;
        s = area_triangle(p2, p3, c) + area_circle_angle(p3, p1, c, r);
        return s * flag;
    }
}
//圆与多边形相交面积
double area_polygon_circle(int n, point p[], CP &c, double r)
{
    double ans = 0;
    for (int i = 0; i < n; ++i)
        ans += area_triangle_circle(c, r, p[i], p[(i + 1) % n]);
    return fabs(ans);
}
//圆与圆相交面积
double area_circle_circle(CC c1, CC c2)
{
    double d = dis(c1.c, c2.c);
    if (c1.r + c2.r <= d)
        return 0;         //两圆相离
    if (c1.r - c2.r >= d) //两园相含 c1包含c2
        return acos(-1.0) * c2.r * c2.r;
    if (c2.r - c1.r >= d) //两园相含 c2包含c1
        return acos(-1.0) * c1.r * c1.r;
    //两圆相交
    double angle1 = acos((c1.r * c1.r + d * d - c2.r * c2.r) / (2 * d * c1.r));
    double angle2 = acos((c2.r * c2.r + d * d - c1.r * c1.r) / (2 * d * c2.r));
    return c1.r * c1.r * angle1 + c2.r * c2.r * angle2 - sin(angle1) * c1.r * d;
}

//球面
//lat表示纬度,-90<=w<=90,lng表示经度
//返回两点所在大圆劣弧对应圆心角,0<=angle<=pi
double angle(double lng1, double lat1, double lng2, double lat2)
{
    double dlng = fabs(lng1 - lng2) * pi / 180;
    while (dlng >= pi + pi)
        dlng -= pi + pi;
    if (dlng > pi)
        dlng = pi + pi - dlng;
    lat1 *= pi / 180;
    lat2 *= pi / 180;
    return acos(cos(lat1) * cos(lat2) * cos(dlng) + sin(lat1) * sin(lat2));
}
//计算距离,r为球半径
double line_dist(double r, double lng1, double lat1, double lng2, double lat2)
{
    double dlng = fabs(lng1 - lng2) * pi / 180;
    while (dlng >= pi + pi)
        dlng -= pi + pi;
    if (dlng > pi)
        dlng = pi + pi - dlng;
    lat1 *= pi / 180;
    lat2 *= pi / 180;
    return r * sqrt(2 - 2 * (cos(lat1) * cos(lat2) * cos(dlng) + sin(lat1) * sin(lat2)));
}
//球面距离
inline double sphere_dist(double r, double lng1, double lat1, double lng2, double lat2)
{
    return r * angle(lng1, lat1, lng2, lat2);
}
//球面三角形余弦定理
//cos c = cos a * cos b + sin a * sin b * cosC

struct gpoint
{
    int x, y;
};
//多边形网格点个数
//一条左开右闭的线段上的格点数为gcd(x2-x1, y2-y1)
int gcd(int a, int b)
{
    return b ? gcd(b, a % b) : a;
}
int grid_onedge(int n, gpoint *p)
{
    int i, ret = 0;
    for (i = 0; i < n; i++)
        ret += gcd(abs(p[i].x - p[(i + 1) % n].x), abs(p[i].y - p[(i + 1) % n].y));
    return ret;
}
int grid_inside(int n, gpoint *p)
{
    int i, ret = 0;
    for (i = 0; i < n; i++)
        ret += p[(i + 1) % n].y * (p[i].x - p[(i + 2) % n].x);
    return (abs(ret) - grid_onedge(n, p)) / 2 + 1;
} //pick 定理

//区域中点集个数
//意义不明

#endif

#ifndef Geometry_3D_H
#define Geometry_3D_H
struct point3
{
    double x, y, z;
    point3()
    {
        x = 0;
        y = 0;
        z = 0;
    }
    point3(double sx, double sy, double sz)
    {
        x = sx;
        y = sy;
        z = sz;
    }
    bool operator<(const point3 &b) const
    {
        if (b.x == x)
        {
            if (y == b.y)
                return z < b.z;
            return y < b.y;
        }
        return x < b.x;
    }
    point3 operator-(const point3 &b)
        const
    {
        point3 a;
        a.x = x - b.x;
        a.y = y - b.y;
        a.z = z - b.z;
        return a;
    }
    point3 operator+(const point3 &b)
        const
    {
        point3 a;
        a.x = x + b.x;
        a.y = y + b.y;
        a.z = z + b.z;
        return a;
    }
    point3 operator/(const double &c)
        const
    {
        point3 a;
        a.x = x / c;
        a.y = y / c;
        a.z = z / c;
        return a;
    }
    point3 operator*(const double &c)
        const
    {
        point3 a;
        a.x = x * c;
        a.y = y * c;
        a.z = z * c;
        return a;
    }
    bool operator==(const point3 &p)
        const
    {
        return zero(x - p.x) && zero(y - p.y) &&
               zero(z - p.z);
    }
};

struct line3
{
    point3 a, b;
    line3(){};
    line3(const point3 &p1, const point3 &p2)
    {
        a = p1;
        b = p2;
    }
};
struct plane3
{
    point3 a, b, c;
};
typedef const point3 CP3;
typedef const line3 CL3;

// 计算cross  product  U x  V
point3 xmult(CP3 &u, CP3 &v)
{
    point3 ret;
    ret.x = u.y * v.z - v.y * u.z;
    ret.y = u.z * v.x - u.x * v.z;
    ret.z = u.x * v.y - u.y * v.x;
    return ret;
}

// 计算dot  product  U . V
double dmult(CP3 &u, CP3 &v)
{
    return u.x * v.x + u.y * v.y + u.z * v.z;
}

// 取平面法向量
point3 pvec(const plane3 &s)
{
    return xmult(s.a - s.b, s.b - s.c);
}
point3 pvec(CP3 &s1, CP3 &s2, CP3 &s3)
{
    return xmult(s1 - s2, s2 - s3);
}

// 两点距离,  单参数取向量大小
double dis(CP3 &p1, CP3 &p2)
{
    return sqrt(sqr(p1.x - p2.x) + sqr(p1.y - p2.y) + sqr(p1.z - p2.z));
}
double dis2(CP3 &p1, CP3 &p2)
{
    return sqr(p1.x - p2.x) + sqr(p1.y - p2.y) + sqr(p1.z - p2.z);
}

// 向量大小
double len(CP3 &p)
{
    return sqrt(sqr(p.x) + sqr(p.y) + sqr(p.z));
}

// 判三点共线
int dots_inline(CP3 &p1, CP3 &p2, CP3 &p3)
{
    return len(xmult(p1 - p2, p2 - p3)) < eps;
}

// 判四点共面
int dots_onplane(CP3 &a, CP3 &b, CP3 &c, CP3 &d)
{
    return zero(dmult(pvec(a, b, c), d - a));
}

// 判点是否在线段上,包括端点和共线
int dot_online_in(CP3 &p, const line3 &l)
{
    return zero(len(xmult(p - l.a, p - l.b))) &&
           (l.a.x - p.x) * (l.b.x - p.x) < eps && (l.a.y - p.y) * (l.b.y - p.y) < eps && (l.a.z - p.z) * (l.b.z - p.z) < eps;
}
int dot_online_in(CP3 &p, CP3 &l1, CP3 &l2)
{
    return zero(len(xmult(p - l1, p - l2))) &&
           (l1.x - p.x) * (l2.x - p.x) < eps && (l1.y - p.y) * (l2.y - p.y) < eps && (l1.z - p.z) * (l2.z - p.z) < eps;
}

// 判点是否在线段上,  不包括端点
int dot_online_ex(CP3 &p, const line3 &l)
{
    return dot_online_in(p, l) && (!(p == l.a)) && (!(p == l.b));
}
int dot_online_ex(CP3 &p, CP3 &l1, CP3 &l2)
{
    return dot_online_in(p, l1, l2) && (!(p == l1)) && (!(p == l2));
}

//判点是否在空间三角形上,包括边界,三点共线无意义
int dot_inplane_in(CP3 &p, const plane3 &s)
{
    return zero(len(xmult(s.a - s.b, s.a - s.c)) - len(xmult(p - s.a, p - s.b)) - len(xmult(p - s.b, p - s.c)) - len(xmult(p - s.c, p - s.a)));
}
int dot_inplane_in(CP3 &p, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return zero(len(xmult(s1 - s2, s1 - s3)) - len(xmult(p - s1, p - s2)) - len(xmult(p - s2, p - s3)) - len(xmult(p - s3, p - s1)));
}

// 判点是否在空间三角形上,不包括边界,三点共线无意义
int dot_inplane_ex(CP3 &p, const plane3 &s)
{
    return dot_inplane_in(p, s) && len(xmult(p - s.a, p - s.b)) > eps &&
           len(xmult(p - s.b, p - s.c)) > eps &&
           len(xmult(p - s.c, p - s.a)) > eps;
}
int dot_inplane_ex(CP3 &p, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return dot_inplane_in(p, s1, s2, s3) &&
           len(xmult(p - s1, p - s2)) > eps &&
           len(xmult(p - s2, p - s3)) > eps &&
           len(xmult(p - s3, p - s1)) > eps;
}

// 判两点在线段同侧,点在线段上返回0,不共面无意义
int same_side(CP3 &p1, CP3 &p2, const line3 &l)
{
    return dmult(xmult(l.a - l.b, p1 - l.b), xmult(l.a - l.b, p2 - l.b)) > eps;
}
int same_side(CP3 &p1, CP3 &p2, CP3 &l1, CP3 &l2)
{
    return dmult(xmult(l1 - l2, p1 - l2), xmult(l1 - l2, p2 - l2)) > eps;
}

//异侧
int opposite_side(CP3 &p1, CP3 &p2, const line3 &l)
{
    return dmult(xmult(l.a - l.b, p1 - l.b), xmult(l.a - l.b, p2 - l.b)) < -eps;
}
int opposite_side(CP3 &p1, CP3 &p2, CP3 &l1, CP3 &l2)
{
    return dmult(xmult(l1 - l2, p1 - l2), xmult(l1 - l2, p2 - l2)) < -eps;
}

// 判两点在平面同侧,  点在平面上返回0
int same_side(CP3 &p1, CP3 &p2, const plane3 &s)
{
    return dmult(pvec(s), p1 - s.a) * dmult(pvec(s), p2 - s.a) > eps;
}
int same_side(CP3 &p1, CP3 &p2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return dmult(pvec(s1, s2, s3), p1 - s1) * dmult(pvec(s1, s2, s3), p2 - s1) > eps;
}

// 判两点在平面异侧,  点在平面上返回0
int opposite_side(CP3 &p1, CP3 &p2, const plane3 &s)
{
    return dmult(pvec(s), p1 - s.a) * dmult(pvec(s), p2 - s.a) < -eps;
}
int opposite_side(CP3 &p1, CP3 &p2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return dmult(pvec(s1, s2, s3), p1 - s1) * dmult(pvec(s1, s2, s3), p2 - s1) < -eps;
}

// 判两直线平行
int parallel(const line3 &u, const line3 &v)
{
    return len(xmult(u.a - u.b, v.a - v.b)) < eps;
}
int parallel(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    return len(xmult(u1 - u2, v1 - v2)) < eps;
}

// 判两平面平行
int parallel(const plane3 &u, const plane3 &v)
{
    return len(xmult(pvec(u), pvec(v))) < eps;
}
int parallel(CP3 &u1, CP3 &u2, CP3 &u3, CP3 &v1, point3 v2, point3 v3)
{
    return len(xmult(pvec(u1, u2, u3), pvec(v1, v2, v3))) < eps;
}

// 判直线与平面平行
int parallel(const line3 &l, const plane3 &s)
{
    return zero(dmult(l.a - l.b, pvec(s)));
}
int parallel(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return zero(dmult(l1 - l2, pvec(s1, s2, s3)));
}

// 判两直线垂直
int perpendicular(const line3 &u, const line3 &v)
{
    return zero(dmult(u.a - u.b, v.a - v.b));
}
int perpendicular(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    return zero(dmult(u1 - u2, v1 - v2));
}

// 判两平面垂直
int perpendicular(const plane3 &u, const plane3 &v)
{
    return zero(dmult(pvec(u), pvec(v)));
}
int perpendicular(CP3 &u1, CP3 &u2, CP3 &u3, CP3 &v1, CP3 &v2, CP3 &v3)
{
    return zero(dmult(pvec(u1, u2, u3), pvec(v1, v2, v3)));
}

//判直线与平面平行
int perpendicular(const line3 &l, const plane3 &s)
{
    return len(xmult(l.a - l.b, pvec(s))) < eps;
}
int perpendicular(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return len(xmult(l1 - l2, pvec(s1, s2, s3))) < eps;
}

// 判两线段相交,包括端点和部分重合
int intersect_in(const line3 &u, const line3 &v)
{
    if (!dots_onplane(u.a, u.b, v.a, v.b))
        return 0;
    if (!dots_inline(u.a, u.b, v.a) || !dots_inline(u.a, u.b, v.b))
        return !same_side(u.a, u.b, v) && !same_side(v.a, v.b, u);
    return dot_online_in(u.a, v) || dot_online_in(u.b, v) || dot_online_in(v.a, u) || dot_online_in(v.b, u);
}
int intersect_in(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    if (!dots_onplane(u1, u2, v1, v2))
        return 0;
    if (!dots_inline(u1, u2, v1) || !dots_inline(u1, u2, v2))
        return !same_side(u1, u2, v1, v2) && !same_side(v1, v2, u1, u2);
    return dot_online_in(u1, v1, v2) ||
           dot_online_in(u2, v1, v2) ||
           dot_online_in(v1, u1, u2) ||
           dot_online_in(v2, u1, u2);
}

// 判两线段相交,不包括端点和部分重合
int intersect_ex(const line3 &u, const line3 &v)
{
    return dots_onplane(u.a, u.b, v.a, v.b) &&
           opposite_side(u.a, u.b, v) &&
           opposite_side(v.a, v.b, u);
}
int intersect_ex(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    return dots_onplane(u1, u2, v1, v2) &&
           opposite_side(u1, u2, v1, v2) &&
           opposite_side(v1, v2, u1, u2);
}

//判线段与空间三角形相交, 包括交于边界和(部分)包含
int intersect_in(const line3 &l, const plane3 &s)
{
    return !same_side(l.a, l.b, s) && !same_side(s.a, s.b, l.a, l.b, s.c) && !same_side(s.b, s.c, l.a, l.b, s.a) && !same_side(s.c, s.a, l.a, l.b, s.b);
}
int intersect_in(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return !same_side(l1, l2, s1, s2, s3) && !same_side(s1, s2, l1, l2, s3) && !same_side(s2, s3, l1, l2, s1) && !same_side(s3, s1, l1, l2, s2);
}

// 判线段与空间三角形相交,  不包括交于边界和(部分) 包含
int intersect_ex(const line3 &l, const plane3 &s)
{
    return opposite_side(l.a, l.b, s) && opposite_side(s.a, s.b, l.a, l.b, s.c) && opposite_side(s.b, s.c, l.a, l.b, s.a) && opposite_side(s.c, s.a, l.a, l.b, s.b);
}
int intersect_ex(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return opposite_side(l1, l2, s1, s2, s3) && opposite_side(s1, s2, l1, l2, s3) && opposite_side(s2, s3, l1, l2, s1) && opposite_side(s3, s1, l1, l2, s2);
}

// 计算两直线交点,注意事先判断直线是否共面和平行
// 线段交点请另外判线段相交同时还是要判断是否平行(!)
point3 intersection(const line3 &u, const line3 &v)
{
    point3 ret = u.a;
    double t = ((u.a.x - v.a.x) * (v.a.y - v.b.y) - (u.a.y - v.a.y) * (v.a.x - v.b.x)) / ((u.a.x - u.b.x) * (v.a.y - v.b.y) - (u.a.y - u.b.y) * (v.a.x - v.b.x));
    ret.x += (u.b.x - u.a.x) * t;
    ret.y += (u.b.y - u.a.y) * t;
    ret.z += (u.b.z - u.a.z) * t;
    return ret;
}

point3 intersection(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    point3 ret = u1;
    double t = ((u1.x - v1.x) * (v1.y - v2.y) - (u1.y - v1.y) * (v1.x - v2.x)) / ((u1.x - u2.x) * (v1.y - v2.y) - (u1.y - u2.y) * (v1.x - v2.x));
    ret.x += (u2.x - u1.x) * t;
    ret.y += (u2.y - u1.y) * t;
    ret.z += (u2.z - u1.z) * t;
    return ret;
}

// 计算直线与平面交点,注意事先判断是否平行,并保证三点不共线!
// 线段和空间三角形交点请另外判断
point3 intersection(const line3 &l, const plane3 &s)
{
    point3 ret = pvec(s);
    double t = (ret.x * (s.a.x - l.a.x) + ret.y * (s.a.y - l.a.y) + ret.z * (s.a.z - l.a.z)) / (ret.x * (l.b.x - l.a.x) + ret.y * (l.b.y - l.a.y) + ret.z * (l.b.z - l.a.z));
    ret.x = l.a.x + (l.b.x - l.a.x) * t;
    ret.y = l.a.y + (l.b.y - l.a.y) * t;
    ret.z = l.a.z + (l.b.z - l.a.z) * t;
    return ret;
}
point3 intersection(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    point3 ret = pvec(s1, s2, s3);
    double t = (ret.x * (s1.x - l1.x) + ret.y * (s1.y - l1.y) + ret.z * (s1.z - l1.z)) /
               (ret.x * (l2.x - l1.x) + ret.y * (l2.y - l1.y) +
                ret.z * (l2.z - l1.z));
    ret.x = l1.x + (l2.x - l1.x) * t;
    ret.y = l1.y + (l2.y - l1.y) * t;
    ret.z = l1.z + (l2.z - l1.z) * t;
    return ret;
}

// 计算两平面交线,  注意事先判断是否平行,  并保证三点不共线 !
line3 intersection(const plane3 &u, const plane3 &v)
{
    line3 ret;
    ret.a = parallel(v.a, v.b, u.a, u.b, u.c) ? intersection(v.b, v.c, u.a, u.b, u.c) : intersection(v.a, v.b, u.a, u.b, u.c);
    ret.b = parallel(v.c, v.a, u.a, u.b, u.c) ? intersection(v.b, v.c, u.a, u.b, u.c) : intersection(v.c, v.a, u.a, u.b, u.c);
    return ret;
}
line3 intersection(CP3 &u1, CP3 &u2, CP3 &u3, CP3 &v1, CP3 &v2, CP3 &v3)
{
    line3 ret;
    ret.a = parallel(v1, v2, u1, u2, u3) ? intersection(v2, v3, u1, u2, u3) : intersection(v1, v2, u1, u2, u3);
    ret.b = parallel(v3, v1, u1, u2, u3) ? intersection(v2, v3, u1, u2, u3) : intersection(v3, v1, u1, u2, u3);
    return ret;
}

// 点到直线距离
double ptoline(CP3 &p, const line3 &l)
{
    return len(xmult(p - l.a, l.b - l.a)) /
           dis(l.a, l.b);
}
double ptoline(CP3 &p, CP3 &l1, CP3 &l2)
{
    return len(xmult(p - l1, l2 - l1)) / dis(l1, l2);
}

// 点到平面距离
double ptoplane(CP3 &p, const plane3 &s)
{
    return fabs(dmult(pvec(s), p - s.a)) /
           len(pvec(s));
}
double ptoplane(CP3 &p, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return fabs(dmult(pvec(s1, s2, s3), p - s1)) / len(pvec(s1, s2, s3));
}

// 直线到直线距离
double linetoline(const line3 &u, const line3 &v)
{
    point3 n = xmult(u.a - u.b, v.a - v.b);
    return fabs(dmult(u.a - v.a, n)) / len(n);
}
double linetoline(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    point3 n = xmult(u1 - u2, v1 - v2);
    return fabs(dmult(u1 - v1, n)) / len(n);
}

// 两直线夹角cos值
double angle_line_line_cos(const line3 &u, const line3 &v)
{
    return dmult(u.a - u.b, v.a - v.b) / len(u.a - u.b) / len(v.a - v.b);
}
double angle_line_line_cos(CP3 &u1, CP3 &u2, CP3 &v1, CP3 &v2)
{
    return dmult(u1 - u2, v1 - v2) / len(u1 - u2) / len(v1 - v2);
}

// 两平面夹角cos值
double angle_plane_plane_cos(const plane3 &u, const plane3 &v)
{
    return dmult(pvec(u), pvec(v)) / len(pvec(u)) / len(pvec(v));
}
double angle_plane_plane_cos(CP3 &u1, CP3 &u2, CP3 &u3, CP3 &v1, CP3 &v2, CP3 &v3)
{
    return dmult(pvec(u1, u2, u3), pvec(v1, v2, v3)) / len(pvec(u1, u2, u3)) / len(pvec(v1, v2, v3));
}

// 直线平面夹角sin值
double angle_line_plane_sin(const line3 &l, const plane3 &s)
{
    return dmult(l.a - l.b, pvec(s)) / len(l.a - l.b) / len(pvec(s));
}
double angle_line_plane_sin(CP3 &l1, CP3 &l2, CP3 &s1, CP3 &s2, CP3 &s3)
{
    return dmult(l1 - l2, pvec(s1, s2, s3)) /
           len(l1 - l2) / len(pvec(s1, s2, s3));
}

//三角形有向面积
double area_triangle(const plane3 &p)
{
    return len(xmult(p.b - p.a, p.c - p.a)) / 2;
}
double area_triangle(CP3 &p1, CP3 &p2, CP3 &p3)
{
    return len(xmult(p2 - p1, p3 - p1)) / 2;
}
double area_triangle(CP3 &p2, CP3 &p3)
{
    return len(xmult(p2, p3)) / 2;
}

//多边形有向面积
double area_polygon(int n, point3 *p)
{
    double s = 0;
    for (int i = 0; i < n; ++i)
        s += len(xmult(p[i], p[(i + 1) % n])) / 2;
    return s;
}

//求四面体有向体积
double volume_tetrahedron(CP3 &p1, CP3 &p2, CP3 &p3, CP3 &p4)
{
    return dmult(xmult(p1 - p4, p2 - p4), p3 - p4) / 6;
}
double volume_tetrahedron(CP3 &p1, CP3 &p2, CP3 &p3)
{
    return dmult(xmult(p1, p2), p3) / 6;
}
double volume_tetrahedron(const plane3 &p)
{
    return dmult(xmult(p.a, p.b), p.c) / 6;
}

//多面体有向体积
double volume_polygon(int n, plane3 *polygon)
{
    double c = 0;
    for (int i = 0; i < n; ++i)
        c += volume_tetrahedron(polygon[i]);
    return c;
}

//三角形重心
point3 barycenter(const point3 &a, const point3 &b, const point3 &c)
{
    return a + b + c / 3;
}

// 四面体重心
point3 barycenter(const point3 &a, const point3 &b, const point3 &c, const point3 &d)
{
    return (a + b) + (c + d) / 4;
}

// 多面体重心
point3 barycenter(int n, plane3 *polygon)
{
    point3 c;
    double v = 0;
    for (int i = 0; i < n; ++i)
    {
        double j = volume_tetrahedron(polygon[i]);
        v += j;
        c = c + (polygon[i].a + polygon[i].b + polygon[i].c) * j;
    }
    return c / (4 * v);
}

//凸包
//不是所有点共面 O(nlgn)
const int MAXN3 = 500;
const int MAXM3 = 250000;
struct NODE
{
    int p[4], next, out;
    point3 f;
} s[MAXM3];
int edge[MAXN3][MAXN3];
int tot;
int next(int x)
{
    if (s[x].next == x)
        return x;
    return s[x].next = next(s[x].next);
}
void add(int a, int b, int c, point3 *p)
{
    s[tot].p[0] = a;
    s[tot].p[1] = b;
    s[tot].p[2] = c;
    s[tot].p[3] = a;
    s[tot].f = xmult(p[b] - p[a], p[c] - p[a]);
    s[tot].out = false;
    for (int i = 0; i < 3; ++i)
        edge[s[tot].p[i]][s[tot].p[i + 1]] = tot;
    ++tot;
}
void add(int a, int b, int c, int d, point3 *p)
{
    point3 f = xmult(p[b] - p[a], p[c] - p[a]);
    if (dmult(f, p[d] - p[a]) > 0)
        add(a, c, b, p);
    else
        add(a, b, c, p);
}
int get_convex(int n, point3 *p, plane3 *convex)
{
    //increment algorithm
    for (int i = 0; i < MAXN; ++i)
        s[i].next = i;
    tot = 0;
    for (int i = 3; i < n; ++i)
        if (!dots_onplane(p[0], p[1], p[2], p[i]))
        {
            std::swap(p[i], p[3]);
            break;
        }
    add(0, 1, 2, 3, p);
    add(2, 3, 0, 1, p);
    add(3, 1, 0, 2, p);
    add(3, 1, 2, 0, p);
    for (int i = 4; i < n; ++i)
    {
        for (int j = next(0); j < tot; j = next(j + 1))
            s[j].out = dmult(s[j].f, p[i] - p[s[j].p[0]]) > 0;
        int c = tot;
        for (int j = next(0); j < tot; j = next(j + 1))
            if (s[j].out)
            {
                for (int k = 0; k < 3; ++k)
                    if (!s[edge[s[j].p[k + 1]][s[j].p[k]]].out)
                        add(s[j].p[k], s[j].p[k + 1], i, p);
                s[j].next = j + 1;
            }
    }
    int i, j;
    for (i = 0, j = next(0); j < tot; ++i, j = next(j + 1))
    {
        convex[i].a = p[s[j].p[0]];
        convex[i].b = p[s[j].p[1]];
        convex[i].c = p[s[j].p[2]];
    }
    return i;
}

#endif
#ifndef Geometry_2D_Complex_H
#define Geometry_2D_Complex_H
typedef std::complex<double> pt;
bool cmp(pt a,pt b){
    return std::make_pair(a.real(),a.imag())<std::make_pair(b.real(),b.imag());
}
#endif
#endif
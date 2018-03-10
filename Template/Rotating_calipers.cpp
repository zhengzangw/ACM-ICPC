#include <iostream>
#include <vector>
#include <cmath>
using namespace std;
const int MAX_N = 1 << 17;
const int eps = 1e-18;
int n,m,t,x,y;
struct point{
    double x,y;
};
point p[MAX_N];

bool cmp(point a, point b){
    if (a.x>b.x) return true;
    else return false;
}

double dist(point a,point b){
    return sqrt((b.x-a.x)*(b.x-a.x)+(b.y-a.y)*(b.y-a.y));
}

double det(double x1,double y1, double x2, double y2){
    return x1*y2-x2*y1;
}

double rotate_calipers(){
    if (n==2) return dist(p[0],p[1]);
    int i=0, j=0;
    for (int k=0;k<n;k++){
        if (!cmp(p[i],p[k])) i=k;
        if (cmp(p[j],p[k])) j=k;
    }
    
    double res=0;
    int si = i, sj = j;
    while (i!=sj||j!=si){
        res = max(res,dist(p[j],p[i]));
        if (det(p[(i+1)%n].x-p[i].x,p[(i+1)%n].y-p[i].y,p[(j+1)%n].x-p[j].x,p[(j+1)%n].y-p[j].y)<0) i = (i+1)%n;
        else j = (j+1)%n;
    }
    
    return res;
}

int main(){
    cin >> n;
    for (int i=0;i<n;i++)//Suppose the input is a convex_hull
        cin >> p[i].x >> p[i].y;
    double ans = rotate_calipers();
    cout << int(ans*ans);
    return 0;
}
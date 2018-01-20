#include <iostream>
#include <vector>
#include <string>
#include <cstdio>
using namespace std;

const int MAX_N = 1000;
int n,m,t;
int heap[MAX_N],f=0;

void push(int x){ //根节点为0
    int i = f++;
    while (i>0){
        int p = (i-1)/2;
        if (heap[p]<=x) break;
        heap[i] = heap[p];
        i = p;
    }
    heap[i] = x;
}

int pop(){
    int ret = heap[0];
    int x = heap[f--];
    int i = 0;
    while(i*2+1 < f){
        int a = i*2+1, b= i*2 + 2;
        if (b < f && heap[b]<heap[a]) a = b;
        if (heap[a] >= x) break;
        heap[i] = heap[a];
        i = a;
    }
    heap[i] = x;
    return ret;
}


int main(){
    cin >> n;
    for (int i=0;i<n;i++){
        cin >> t;
        push(t);
    }
    for (int i=0;i<n;i++){
        cout << pop();
    }
    return 0;
}
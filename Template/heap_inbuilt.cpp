#include <iostream>
#include <queue>
using namespace std;

const int MAX_N = 1000;
int n,m,t;
int heap[MAX_N],f=0;


int main(){
    cin >> n;
    priority_queue<int, vector<int>, less<int>> pque;
    // less 大根堆 greater 小根堆
    for (int i=0;i<n;i++){
        cin >> t;
        pque.push(t);
    }
    while (!pque.empty()){
        cout << pque.top();
        pque.pop();
    }
    return 0;
}
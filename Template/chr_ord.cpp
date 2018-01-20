#include <iostream>
using namespace std;
const int MAX_N = 1 << 17;
const int INF = 1 << 27;
int n,m,t;

int main(){
    char c;
    cin >> c;
    cout << (char)(((int)c - (int)'a' + 2)%26 + (int)'a' - 1);
    return 0;
}
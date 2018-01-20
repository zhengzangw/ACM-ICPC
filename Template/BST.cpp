#include <iostream>
using namespace std;
int n,m,t;

struct node{
	int val;
	node *lch, *rch;
};

node *insert(node *p,int x){
	if (p == NULL) {
        node *q = new node;
        q->val = x;
        q->lch = q->rch = NULL;
        return q;
    }
    else {
        if (x < p->val) p->lch = insert(p->lch,x);
        else p->rch = insert(p->rch,x);
        return p;
    }
}

bool find(node *p, int x){
    if (p == NULL) return false;
    if (x == p->val) return true;
    if (x < p->val) return find(p->lch,x);
     else return find(p->rch,x);
}

node *remove(node *p, int x){
    if (p == NULL) return NULL;
    else if (x < p->val) p->lch = remove(p->lch,x);
    else if (x > p->val) p->rch = remove(p->rch,x);

    else if (p->lch == NULL) {
        node *q = p->rch;
        delete p;
        return q;
    }

    else if (p->lch->rch == NULL){
        node *q = p->lch;
        q->rch = p->rch;
        delete p;
        return q;
    }
    else {
        node *q;
        for (q = p->lch; q->rch->rch != NULL; q = q->rch);
        node *r = q->rch;
        q->rch = r->lch;
        r->lch = p->lch;
        r->rch = p->rch;
        delete p;
        delete r;
    }
    return p;
}

int main(){
    cin >> n;

    node *root = NULL;
    

    for (int i=0;i<n;i++){
        cin >> t;
        root = insert(root,t);
    }
    
    cin >> t;
    char c;
    while (t--){
        cin >> c >> m;
        if (c == 'd') root = remove(root,m);
        if (c == 'q') { 
            if (find(root,m)) cout << "IN";
            else cout << "NO!";
        }
    }
    return 0;
}
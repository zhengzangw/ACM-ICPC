#include <iostream>
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <string>
#include <cstring>
#include <map>
#include <set>
#include <vector>
#include <queue>
#include <stack>
#include <iomanip>
#include <cstdlib>
#include <climits>
//#include <unordered_map>
typedef long long LL;
typedef unsigned int UI;
typedef unsigned long long ULL;
const LL mod = 1000000007;
#define MAXN 100
#define eps 1e-10
#define fabs ((x) > 0 ? (x) : -(x))
#define random(x) (rand() % x)
using namespace std;

class BST{
private:
    struct node{
        int key;
        node *lch,*rch,*p;
        node(){
            lch = NULL; rch = NULL; p = NULL;
        }
    };
    int firstprint;
    node *root;
    node* SEARCH(node *x,int k){
        while (x!=NULL&&k!=x->key){
            if (k<x->key) x=x->lch;
            else x=x->rch;
        }
        return x;
    }
    node* MINIMUM(node *x){
        while (x->lch!=NULL)
            x = x->lch;
        return x;
    }
    node* MAXIMUN(node *x){
        while (x->rch!=NULL)
            x = x->rch;
        return x;
    }
    node* SUCCESSOR(node *x){
        if (x->rch!=NULL)
            return MINIMUM(x->rch);
        node* y = x->p;
        while (y!=NULL&&x==y->rch){
            x = y;
            y = y->p;
        }
        return y;
    }
    node* PRECCESSOR(node *x){
        if (x->lch!=NULL)
            return MAXIMUN(x->lch);
        node* y = x->p;
        while (y!=NULL&&x==y->lch){
            x = y;
            y = y->p;
        }
        return y;
    }
    void INSERT(node*z){
        node *y = NULL;
        node *x = root;
        while (x!=NULL){
            y = x;
            if (z->key<x->key) x = x->lch;
            else x = x->rch;
        }
        z->p = y;
        if (y==NULL) root = z;
        else if (z->key<y->key) y->lch = z;
        else y->rch = z;
    }
    void TRANSPLANT(node *u,node *v){
        if (u->p==NULL)
            root = v;
        else if (u==u->p->lch)
            u->p->lch = v;
        else u->p->rch = v;
        if (v!=NULL)
            v->p = u->p;
    }
    void DELETE(node *z){
        if (z->lch==NULL)
            TRANSPLANT(z,z->rch);
        else if (z->rch==NULL)
            TRANSPLANT(z,z->lch);
        else{
            node *y = MINIMUM(z->rch);
            if (y->p!=z){
                TRANSPLANT(y,y->rch);
                y->rch = z->rch;
                y->rch->p = y;
            }
            TRANSPLANT(z,y);
            y->lch = z->lch;
            y->lch->p = y;
        }
    }
    void PREWALK(node* x){
        if (x!=NULL){
            if (firstprint){
                firstprint = 0;
                cout << x->key;
            } else cout << ' '<< x->key;
            PREWALK(x->lch);
            PREWALK(x->rch);
        }
    }
    void INWALK(node* x){
        if (x!=NULL){
            INWALK(x->lch);
            if (firstprint)
            {
                firstprint = 0;
                cout << x->key;
            }
            else
                cout << ' ' << x->key;
            INWALK(x->rch);
        }
        
    }
    void POSTWALK(node* x){
        if (x!=NULL){
            POSTWALK(x->lch);
            POSTWALK(x->rch);
            if (firstprint)
            {
                firstprint = 0;
                cout << x->key;
            }
            else
                cout << ' ' << x->key;
        }
    }
    
public:
    BST(){
        root = NULL;
    }
    node *TREE_ROOT(){
        return root;
    }
    int TREE_SEARCH(int k){
        node* temp = SEARCH(root,k);
        if (temp==NULL) return 0;
        else return 1;
    }
    void TREE_INSEART(int k){
        node* temp = new node;
        temp -> key = k;
        INSERT(temp);
    }
    void TREE_DELETE(int k){
        node* temp = SEARCH(root,k);
        DELETE(temp);
    }
    void Preord_walk(){
        firstprint = 1;
        PREWALK(root);
    }
    void Inord_walk(){
        firstprint = 1;
        INWALK(root);
    }
    void Postord_walk(){
        firstprint = 1;
        POSTWALK(root);
    }
};

int main()
{
    int n,k;
    char op;
    cin >> n;
    BST bst;
    for (int i=0;i<n;i++){
        cin >> op >> k;
        switch(op){
            case 'A': bst.TREE_INSEART(k);
            break;
            case 'D': bst.TREE_DELETE(k);
            break;
        } 
    }
    bst.Inord_walk(); cout << endl;
    bst.Preord_walk(); cout << endl;
    bst.Postord_walk();
    return 0;
}
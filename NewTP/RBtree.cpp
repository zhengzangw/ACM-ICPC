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

class RB_BST
{
  private:
    struct node
    {
        bool NIL;
        int key;
        char color;
        node *lch, *rch, *p;
        node(int val)
        {
            key = val;
            NIL = false;
            lch = NULL;
            rch = NULL;
            p = NULL;
            color = 'R';
        }
        node()
        {
            NIL = true;
            lch = NULL;
            rch = NULL;
            p = NULL;
            color = 'B';
        }
    };
    int firstprint;
    node *root;
    node *NIL;
    void LEFT_ROTATION(node *x)
    {
        node *y = x->rch;
        x->rch = y->lch;
        if (y->lch != NIL)
            y->lch->p = x;
        y->p = x->p;
        if (x->p == NIL)
            root = y;
        else if (x == x->p->lch)
            x->p->lch = y;
        else
            x->p->rch = y;
        y->lch = x;
        x->p = y;
    }
    void RIGHT_ROTATION(node *x)
    {
        node *y = x->lch;
        x->lch = y->rch;
        if (y->rch != NIL)
            y->rch->p = x;
        y->p = x->p;
        if (x->p == NIL)
            root = y;
        else if (x == x->p->lch)
            x->p->lch = y;
        else
            x->p->rch = y;
        y->rch = x;
        x->p = y;
    }
    node *SEARCH(node *x, int k)
    {
        while (x != NIL && k != x->key)
        {
            if (k < x->key)
                x = x->lch;
            else
                x = x->rch;
        }
        return x;
    }
    int SEARCH_count(node *x, int k)
    {
        int cnt = 0;
        while (x != NIL && k != x->key)
        {
            cnt++;
            if (k < x->key)
                x = x->lch;
            else
                x = x->rch;
        }
        return cnt;
    }
    node *MINIMUM(node *x)
    {
        while (x->lch != NIL)
            x = x->lch;
        return x;
    }
    node *MAXIMUN(node *x)
    {
        while (x->rch != NIL)
            x = x->rch;
        return x;
    }
    node *SUCCESSOR(node *x)
    {
        if (x->rch != NIL)
            return MINIMUM(x->rch);
        node *y = x->p;
        while (y != NIL && x == y->rch)
        {
            x = y;
            y = y->p;
        }
        return y;
    }
    node *PRECCESSOR(node *x)
    {
        if (x->lch != NIL)
            return MAXIMUN(x->lch);
        node *y = x->p;
        while (y != NIL && x == y->lch)
        {
            x = y;
            y = y->p;
        }
        return y;
    }
    void TRANSPLANT(node *u, node *v)
    {
        if (u->p == NIL)
            root = v;
        else if (u == u->p->lch)
            u->p->lch = v;
        else
            u->p->rch = v;
        v->p = u->p;
    }
    void INSERT_FIXUP(node *z)
    {
        while (z->p->color == 'R')
            if (z->p == z->p->p->lch)
            {
                node *y = z->p->p->rch;
                if (y->color == 'R')
                {
                    z->p->color = 'B';
                    y->color = 'B';
                    z->p->p->color = 'R';
                    z = z->p->p;
                }
                else
                {
                    if (z == z->p->rch)
                    {
                        z = z->p;
                        LEFT_ROTATION(z);
                    }
                    z->p->color = 'B';
                    z->p->p->color = 'R';
                    RIGHT_ROTATION(z->p->p);
                }
            }
            else
            {
                node *y = z->p->p->lch;
                if (y->color == 'R')
                {
                    z->p->color = 'B';
                    y->color = 'B';
                    z->p->p->color = 'R';
                    z = z->p->p;
                }
                else
                {
                    if (z == z->p->lch)
                    {
                        z = z->p;
                        RIGHT_ROTATION(z);
                    }
                    z->p->color = 'B';
                    z->p->p->color = 'R';
                    LEFT_ROTATION(z->p->p);
                }
            }
        root->color = 'B';
    }
    void INSERT(node *z)
    {
        node *y = NIL;
        node *x = root;
        while (x != NIL)
        {
            y = x;
            if (z->key < x->key)
                x = x->lch;
            else
                x = x->rch;
        }
        z->p = y;
        if (y == NIL)
            root = z;
        else if (z->key < y->key)
            y->lch = z;
        else
            y->rch = z;
        z->lch = NIL;
        z->rch = NIL;
        z->color = 'R';
        INSERT_FIXUP(z);
    }
    void DELETE_FIXUP(node *x)
    {
        while (x != root && x->color == 'B')
        {
            if (x == x->p->lch)
            {
                node *w = x->p->rch;
                if (w->color == 'R')
                {
                    w->color = 'B';
                    x->p->color = 'R';
                    LEFT_ROTATION(x->p);
                    w = x->p->rch;
                }
                if (w->lch->color == 'B' && w->rch->color == 'B')
                {
                    w->color = 'R';
                    x = x->p;
                }
                else
                {
                    if (w->rch->color == 'B')
                    {
                        w->lch->color = 'B';
                        w->color = 'R';
                        RIGHT_ROTATION(w);
                        w = x->p->rch;
                    }
                    w->color = x->p->color;
                    x->p->color = 'B';
                    w->rch->color = 'B';
                    LEFT_ROTATION(x->p);
                    x = root;
                }
            }
            else
            {
                node *w = x->p->lch;
                if (w->color == 'R')
                {
                    w->color = 'B';
                    x->p->color = 'R';
                    RIGHT_ROTATION(x->p);
                    w = x->p->lch;
                }
                if (w->rch->color == 'B' && w->lch->color == 'B')
                {
                    w->color = 'R';
                    x = x->p;
                }
                else
                {
                    if (w->lch->color == 'B')
                    {
                        w->rch->color = 'B';
                        w->color = 'R';
                        LEFT_ROTATION(w);
                        w = x->p->lch;
                    }
                    w->color = x->p->color;
                    x->p->color = 'B';
                    w->lch->color = 'B';
                    RIGHT_ROTATION(x->p);
                    x = root;
                }
            }
        }
        x->color = 'B';
    }

    void DELETE(node *z)
    {
        node *y = z;
        node *x;
        char y_original_color = y->color;
        if (z->lch == NIL)
        {
            x = z->rch;
            TRANSPLANT(z, z->rch);
        }
        else if (z->rch == NIL)
        {
            x = z->lch;
            TRANSPLANT(z, z->lch);
        }
        else
        {
            y = MINIMUM(z->rch);
            y_original_color = y->color;
            x = y->rch;
            if (y->p == z)
                x->p = y;
            else
            {
                TRANSPLANT(y, y->rch);
                y->rch = z->rch;
                y->rch->p = y;
            }
            TRANSPLANT(z, y);
            y->lch = z->lch;
            y->lch->p = y;
            y->color = z->color;
        }
        if (y_original_color == 'B')
            DELETE_FIXUP(x);
    }
    void PREWALK(node *x)
    {
        if (x!=NIL)
        {
            cout << x->key << ':';
            if (x->color == 'B')
                cout << "black\n";
            else
                cout << "red\n";
            PREWALK(x->lch);
            PREWALK(x->rch);
        }
    }
    void INWALK(node *x)
    {
        if (x!=NIL)
        {
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
    void POSTWALK(node *x)
    {
        if (x!=NIL)
        {
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
    RB_BST()
    {
        NIL = new node();
        root = NIL;
    }
    node *TREE_ROOT()
    {
        return root;
    }
    int TREE_SEARCH_count(int k)
    {
        int cnt = SEARCH_count(root, k);
        return cnt;
    }
    void TREE_INSEART(int k)
    {
        node *temp = new node(k);
        temp->lch = NIL; temp->rch = NIL;
        temp->p = NIL;
        INSERT(temp);
    }
    void TREE_DELETE(int k)
    {
        node *temp = SEARCH(root, k);
        DELETE(temp);
    }
    void Preord_walk()
    {
        firstprint = 1;
        PREWALK(root);
    }
    void Inord_walk()
    {
        firstprint = 1;
        INWALK(root);
    }
    void Postord_walk()
    {
        firstprint = 1;
        POSTWALK(root);
    }
};

int main()
{
    int n, k, m, cnt;
    cin >> n >> m;
    RB_BST bst;
    for (int i = 0; i < n; i++)
    {
        cin >> k;
        bst.TREE_INSEART(k);
    }
    bst.Preord_walk();
    int firstprint = 1;
    for (int i = 0; i < m; i++)
    {
        cin >> k;
        cnt = bst.TREE_SEARCH_count(k);
        if (firstprint)
        {
            cout << cnt;
            firstprint = 0;
        }
        else
            cout << endl
                 << cnt;
    }
    return 0;
}

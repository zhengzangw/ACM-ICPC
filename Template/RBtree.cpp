#ifndef RBBST_H
#define RBBST_H

#include <cmath>
#include <cstdio>
template <typename T>
class RBT
{
    protected:
        struct RBTnode;
        int len;
        int llen; //repeat
        RBTnode *_root;
        RBTnode *_hot;
        RBTnode *_NIL;
        void left_rotate(RBTnode *);
        void right_rotate(RBTnode *);
        void transplant(RBTnode *, RBTnode *);
        RBTnode *find(T, const int);
        RBTnode *findkth(int, RBTnode *);
        int find_rank(T, RBTnode *);

        void insert_fixup(RBTnode*);
        void delete_fixup(RBTnode*);

    public:
        struct iterator;
        RBT(): _root(NULL), len(0), llen(0) {}

        iterator insert(T);
        bool remove(T);

        int get_rank(T);
        iterator kth(int);
        iterator lower_bound(T);
        iterator upper_bound(T);
        iterator find(T);
        int count(T);

        bool empty();
        int size();
};

template <typename T>
struct RBT<T>::RBTnode
{
    bool NIL;
    int value;
    char color;
    RBTnode *lch, *rch, *p;
    RBTnode(int val) : value(val),NIL(false),lch(NULL),rch(NULL),p(NULL),color('R') {}
    RBTnode(): NIL(true),lch(NULL),rch(NULL),p(NULL) {}
    RBTnode *succ(){
        RBTnode *x = this;
        if (x->rch != NULL)
        {
            while (x->lch != NULL)
                x = x->lch;
            return x;
        }
        while (x->p != NULL && x == x->p->rch)
            x = x->p;
        return x->p;
    }
    RBTnode *pre()
    {
        RBTnode *x = this;
        if (x->lch != NULL)
        {
            while (x->rch != NULL)
                x = x->rch;
            return x;
        }
        while (x->p != NULL && x == x->p->lch)
            x = x->p;
        return x->p;
    }
};

template <typename T>
struct RBT<T>::iterator
{
  private:
    RBTnode *_real_node;

  public:
    iterator &operator++()
    {
        _real_node = _real_node->succ();
        return *this;
    }
    iterator &operator--()
    {
        _real_node = _real_node->pre();
        return *this;
    }
    T operator*()
    {
        assert(_real_node != NULL);
        return _real_node->value;
    }
    iterator(RBTnode *node_nn = NULL) : _real_node(node_nn) {}
    iterator(iterator const &iter) : _real_node(iter._real_node) {}
};

template <typename T>
void RBT<T>::left_rotate(RBTnode *x)
{
    RBTnode *y = x->rch;
    x->rch = y->lch;
    if (y->lch != _NIL)
        y->lch->p = x;
    y->p = x->p;
    if (x->p == _NIL)
        _root = y;
    else if (x == x->p->lch)
        x->p->lch = y;
    else
        x->p->rch = y;
    y->lch = x;
    x->p = y;
}

template <typename T>
void RBT<T>::right_rotate(RBTnode *x)
{
    RBTnode *y = x->lch;
    x->lch = y->rch;
    if (y->rch != _NIL)
        y->rch->p = x;
    y->p = x->p;
    if (x->p == _NIL)
        _root = y;
    else if (x == x->p->lch)
        x->p->lch = y;
    else
        x->p->rch = y;
    y->rch = x;
    x->p = y;
}

template <typename T>
typename RBT<T>::RBTnode *RBT<T>::find(T v, const int op)
{
    RBTnode *ptn = _root;
    _hot = NULL;
    while (ptn != NULL && ptn->value != v)
    {
        _hot = ptn;
        ptn->s += op;
        if (ptn->value > v)
            ptn = ptn->lch;
        else
            ptn = ptn->rch;
    }
    return ptn;
}

template <typename T>
void RBT<T>::transplant(RBTnode *u, RBTnode *v)
{
    if (u->p == _NIL)
        _root = v;
    else if (u == u->p->lch)
        u->p->lch = v;
    else
        u->p->rch = v;
    v->p = u->p;
}

template <typename T>
void RBT<T>::insert_fixup(RBTnode *z)
{
    while (z->p->color == 'R')
        if (z->p == z->p->p->lch)
        {
            RBTnode *y = z->p->p->rch;
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
                    left_rotate(z);
                }
                z->p->color = 'B';
                z->p->p->color = 'R';
                right_rotate(z->p->p);
            }
        }
        else
        {
            RBTnode *y = z->p->p->lch;
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
                    right_rotate(z);
                }
                z->p->color = 'B';
                z->p->p->color = 'R';
                left_rotate(z->p->p);
            }
        }
    _root->color = 'B';
}

template <typename T>
typename
RBT<T>::iterator RBT<T>::insert(T k)
{
    llen++;
    RBTnode *tmp = new RBTnode(k);
    if (_root == _NIL)
    {
        _root = tmp;
        len++;
        return iterator(_root);
    }
    RBTnode *ptn = find(k, 1);
    if (ptn != _NIL)
    {
        ptn->count++;
        return iterator(ptn);
    }
    len++;
    if (_hot->value <= k)
        _hot->rch = tmp;
    else
        _hot->lch = tmp;
    tmp->p = _hot;
    insert_fixup(tmp);
}

template <typename T>
void RBT<T>::delete_fixup(RBTnode *x)
{
    while (x != _root && x->color == 'B')
    {
        if (x == x->p->lch)
        {
            RBTnode *w = x->p->rch;
            if (w->color == 'R')
            {
                w->color = 'B';
                x->p->color = 'R';
                left_rotate(x->p);
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
                    right_rotate(w);
                    w = x->p->rch;
                }
                w->color = x->p->color;
                x->p->color = 'B';
                w->rch->color = 'B';
                left_rotate(x->p);
                x = root;
            }
        }
        else
        {
            RBTnode *w = x->p->lch;
            if (w->color == 'R')
            {
                w->color = 'B';
                x->p->color = 'R';
                right_rotate(x->p);
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
                    left_rotate(w);
                    w = x->p->lch;
                }
                w->color = x->p->color;
                x->p->color = 'B';
                w->lch->color = 'B';
                right_rotate(x->p);
                x = _root;
            }
        }
    }
    x->color = 'B';
}

template <typename T>
bool RBT<T>::remove(T k)
{
    RBTnode *z = find(k, -1);
    if (!z)
    {
        find(k, 1);
        //printf("No finded!\n");
        return false;
    }
    llen--;
    if (z->count > 1)
    {
        find(k, 1);
        z->count--;
        return true;
    }
    len--;
    RBTnode *y = z;
    RBTnode *x;
    char y_original_color = y->color;
    if (z->lch == _NIL)
    {
        x = z->rch;
        TRANSPLANT(z, z->rch);
    }
    else if (z->rch == _NIL)
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

template <typename T>
int RBT<T>::size()
{
    return len;
}
template <typename T>
bool RBT<T>::empty()
{
    return _root == NULL;
}

template <typename T>
typename RBT<T>::iterator RBT<T>::lower_bound(T v)
{
    RBTnode *ptn = _root;
    _hot = NULL;
    while (ptn)
    {
        _hot = ptn;
        if (ptn->value < v)
            ptn = ptn->rch;
        else
            ptn = ptn->lch;
    }
    if (_hot->value < v)
        ptn = _hot;
    else
        ptn = _hot->pre();
    return iterator(ptn);
}

template <typename T>
typename RBT<T>::iterator RBT<T>::upper_bound(T v)
{
    RBTnode *ptn = _root;
    _hot = NULL;
    while (ptn)
    {
        _hot = ptn;
        if (ptn->value > v)
            ptn = ptn->lch;
        else
            ptn = ptn->rch;
    }
    if (_hot->value > v)
        ptn = _hot;
    else
        ptn = _hot->succ();
    return iterator(ptn);
}

template <typename T>
typename RBT<T>::iterator RBT<T>::kth(int rank)
{
    assert(rank <= len);
    return iterator(findkth(rank, _root));
}

template <typename T>
typename RBT<T>::RBTnode *RBT<T>::findkth(int rank, RBTnode *ptn)
{
    if (!(ptn->lch))
    {
        if (rank == 1)
        {
            return ptn;
        }
        else
        {
            return findkth(rank - 1, ptn->rch);
        }
    }
    else
    {
        if (ptn->lch->s == rank - 1)
            return ptn;
        else if (ptn->lch->s >= rank)
            return findkth(rank, ptn->lch);
        else
            return findkth(rank - (ptn->lch->s) - 1, ptn->rch);
    }
}

template <typename T>
int RBT<T>::get_rank(T v)
{
    return find_rank(v, _root);
}

template <typename T>
int RBT<T>::find_rank(T v, RBTnode *ptn)
{
    if (!ptn)
        return 1;
    else if (ptn->value >= v)
        return find_rank(v, ptn->lch);
    else
        return (1 + ((ptn->lch) ? (ptn->lch->s) : 0) + find_rank(v, ptn->rch));
}

template <typename T>
int RBT<T>::count(T v)
{
    return find(v, 0)->count;
}

#endif
#ifndef BST_H
#define BST_H
#include <cmath>
#include <cstdio>
#include <cassert>

template <typename T>
class BST
{
  protected:
    struct BSTnode;
    int len;
    int llen; //repeat
    BSTnode *_hot;
    BSTnode *_root;
    void transplant(BSTnode *, BSTnode *);
    BSTnode *find(T, const int);
    BSTnode *findkth(int, BSTnode *);
    int find_rank(T, BSTnode *);

  public:
    struct iterator;
    BST() : _root(NULL), len(0), llen(0) {}

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
struct BST<T>::BSTnode
{
    T value;
    BSTnode *lch, *rch, *p;
    int s, count;
    BSTnode() : lch(NULL), rch(NULL), p(NULL), count(1), s(1) {}
    BSTnode(int x) : value(x), lch(NULL), rch(NULL), p(NULL), count(1), s(1) {}
    BSTnode *succ()
    {
        BSTnode *x = this;
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
    BSTnode *pre()
    {
        BSTnode *x = this;
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
typename BST<T>::BSTnode *BST<T>::find(T v, const int op)
{
    BSTnode *ptn = _root;
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
typename BST<T>::iterator BST<T>::insert(T k)
{
    llen++;
    BSTnode *tmp = new BSTnode(k);
    if (_root == NULL)
    {
        _root = tmp;
        len++;
        return iterator(_root);
    }
    BSTnode *ptn = find(k, 1);
    if (ptn != NULL)
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
    return iterator(tmp);
}

template <typename T>
void BST<T>::transplant(BSTnode *u, BSTnode *v)
{
    if (u->p == NULL)
        _root = v;
    else if (u == u->p->lch)
        u->p->lch = v;
    else
        u->p->rch = v;
    if (v != NULL)
        v->p = u->p;
}

template <typename T>
bool BST<T>::remove(T k)
{
    BSTnode *z = find(k, -1);
    if (!z)
    {
        find(k, 1);
        printf("No finded!\n");
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
    if (z->lch == NULL)
        transplant(z, z->rch);
    else if (z->rch == NULL)
        transplant(z, z->lch);
    else
    {
        BSTnode *y = z->rch;
        while (y->rch != NULL)
            y = y->rch;
        if (y->p != z)
        {
            transplant(y, y->rch);
            y->rch = z->rch;
            y->rch->p = y;
        }
        transplant(z, y);
        y->lch = z->lch;
        y->lch->p = y;
    }
    return true;
}

template <typename T>
int BST<T>::size()
{
    return len;
}
template <typename T>
bool BST<T>::empty()
{
    return _root == NULL;
}

template <typename T>
typename BST<T>::iterator BST<T>::lower_bound(T v)
{
    BSTnode *ptn = _root;
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
typename BST<T>::iterator BST<T>::upper_bound(T v)
{
    BSTnode *ptn = _root;
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
typename BST<T>::iterator BST<T>::kth(int rank)
{
    assert(rank <= len);
    return iterator(findkth(rank, _root));
}

template <typename T>
typename BST<T>::BSTnode *BST<T>::findkth(int rank, BSTnode *ptn)
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
int BST<T>::get_rank(T v)
{
    return find_rank(v, _root);
}

template <typename T>
int BST<T>::find_rank(T v, BSTnode *ptn)
{
    if (!ptn)
        return 1;
    else if (ptn->value >= v)
        return find_rank(v, ptn->lch);
    else
        return (1 + ((ptn->lch) ? (ptn->lch->s) : 0) + find_rank(v, ptn->rch));
}

template <typename T>
int BST<T>::count(T v)
{
    return find(v,0)->count;
}

template <typename T>
struct BST<T>::iterator
{
  private:
    BSTnode *_real_node;

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
    iterator(BSTnode *node_nn = NULL) : _real_node(node_nn) {}
    iterator(iterator const &iter) : _real_node(iter._real_node) {}
};

template <typename T>
typename BST<T>::iterator BST<T>::find(T k)
{
    return iterator(rfind(k, 0));
}
#endif
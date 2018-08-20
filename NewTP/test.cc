#include "BigNum.h";
#include "header.h";
#include "BST.h"

int main()
{
    int n, k;
    char op;
    cin >> n;
    BST bst;
    for (int i = 0; i < n; i++)
    {
        cin >> op >> k;
        switch (op)
        {
        case 'A':
            bst.insert(k);
            break;
        case 'D':
            bst.del(k);
            break;
        }
    }
    bst.inwalk(bst.root);
    cout << endl;
    bst.prewalk(bst.root);
    cout << endl;
    bst.postwalk(bst.root);
    return 0;
}
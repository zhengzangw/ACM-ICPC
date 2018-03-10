#include <iostream>
#include <cstdio>
using namespace std;
#define MAX 26
int n,m,t;

typedef struct TrieNode{
    bool isStr;
    struct TrieNode *next[MAX];
} Trie;

void insert(Trie *root, const char* s)
{
    if(root==NULL||*s=='\0') return;
    Trie *p=root;
    while (*s!='\0')
    {
        if (p->next[*s-'a']==NULL)
        {
            Trie *temp=(Trie *)malloc(sizeof(Trie));
            for (int i = 0;i<MAX;i++) temp->next[i]=NULL;
            temp->isStr=false;
            p->next[*s-'a'] = temp;
            p=p->next[*s-'a'];
        }
        else	p=p->next[*s-'a'];
        s++;
    }
    p -> isStr = true;
}

bool search(Trie *root,const char *s)
{
    Trie *p=root;
    while (p!=NULL&&*s!='\0')
    {
        p=p->next[*s-'a'];
        s++;
    }
    return (p!=NULL&&p->isStr==true);
}

void del(Trie *root)
{
    for (int i=0;i<MAX;i++){
        if(root->next[i]!=NULL) del(root->next[i]);
    }
    free(root);
}


int main(int argc, char *argv[])
{
    int n,m;
    char s[100];
    Trie *root = (Trie *)malloc(sizeof(Trie));
    for (int i=0;i<MAX;i++) root->next[i]=NULL;
    root->isStr=false;
    scanf("%d",&n);
    for (int i=0;i<n;i++){
        scanf("%s",s);
        insert(root,s);
    }
    scanf("%d",&m);
    for (int i=0;i<m;i++){
        scanf("%s",s);
        if (search(root,s)) printf("YES\n"); else printf("NO\n");
    }
    del(root);
    return 0;
}
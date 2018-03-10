#include <iostream>
#include <cstdio>
#include <queue>
using namespace std;
#define MAX 26
int n,m,t,cnt;

struct node{
    int sum;
    node *nxt[MAX];
    node *fail;
};
node* root;

void insert(node*p,const char* s)
{
    for (int i=0;s[i];i++){
        int x=s[i]-'a';
        if (p->nxt[x]==NULL){
            node *newnode=(struct node*)malloc(sizeof(node));
            for (int j=0;j<MAX;j++) newnode->nxt[j]=NULL;
            newnode->sum=0; newnode->fail=NULL;
            p->nxt[x]=newnode;
        }
        p = p->nxt[x];
    }
    p->sum++;
}

void fail_pointer()
{
    queue<node*> q;
    while (!q.empty()) q.pop();
    q.push(root);
    node *p,*temp;
    while (!q.empty()){
        temp=q.front(); q.pop();
        for (int i=0;i<MAX;i++)
            if (temp->nxt[i])
            {
                if (temp == root) temp->nxt[i]->fail=root;
                else {
                    p = temp->fail;
                    while (p){
                        if (p->nxt[i]){
                            temp->nxt[i]->fail=p->nxt[i];
                            break;
                        }
                        p=p->fail;
                    }
                    if (p==NULL) temp->nxt[i]->fail = root;
                }
                q.push(temp->nxt[i]);
            }
    }
}

void AC(node* p,char *ch)
{
    int len=(int)strlen(ch);
    for (int i=0;i<len;i++){
        int x=ch[i]-'a';
        while (!p->nxt[x] && p!= root) p = p->fail;
        p=p->nxt[x];
        if (!p) p=root;
        node* temp=p;
        while (temp!=root){
            if(temp->sum >=0){
                cnt += temp->sum;
                temp->sum=-1;
            } else break;
            temp = temp->fail;
        }
    }
}

int main(int argc, char *argv[])
{
    char s[100];
    cin >> n;
    
    root=(struct node*)malloc(sizeof(node));
    for (int j=0;j<MAX;j++) root->nxt[j]=NULL;
    root->sum=0; root->fail=NULL;
    
    for (int i=0;i<n;i++){
        scanf("%s",s);
        insert(root,s);
    }
    fail_pointer();
    scanf("%s",s);
    cnt = 0;
    AC(root,s);
    cout << cnt;
    return 0;
}

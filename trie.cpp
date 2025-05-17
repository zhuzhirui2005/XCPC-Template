struct trie{
    int siz;
    trie *son[2];
    inline trie(){siz=0,son[0]=son[1]=NULL;}
};
void insert(int dep,trie *p,int k){
    ++p->siz;
    if(dep<0)return;
    int nxt=k>>dep&1;
    if(!p->son[nxt])p->son[nxt]=new trie();
    insert(dep-1,p->son[nxt],k);
}
int query(int dep,trie *p,int k,int lim){
    if(!p)return 0;
    if(dep<0)return p->siz;
    int nxt=k>>dep&1;
    if(lim>>dep&1)return (p->son[nxt]?p->son[nxt]->siz:0)+query(dep-1,p->son[nxt^1],k,lim);
    return query(dep-1,p->son[nxt],k,lim);
}
void insert(int dep,trie *p1,trie *p2,int k){
	if(p1)p2->siz=p1->siz;
	++p2->siz;
	if(dep<0)return;
	int nxt=k>>dep&1;
	if(p1)p2->son[nxt^1]=p1->son[nxt^1];
	p2->son[nxt]=new trie();
	insert(dep-1,p1?p1->son[nxt]:NULL,p2->son[nxt],k);
}
int query(int dep,trie *p1,trie *p2,int k){
	if(dep<0)return 0;
	int nxt=k>>dep&1;
	if(p2->son[nxt^1]&&(!p1||!p1->son[nxt^1]||p2->son[nxt^1]->siz>p1->son[nxt^1]->siz))return query(dep-1,p1?p1->son[nxt^1]:NULL,p2->son[nxt^1],k)|(1<<dep);
	return query(dep-1,p1?p1->son[nxt]:NULL,p2->son[nxt],k);
}
template<class T>
inline V<pii> cart_seq(const V<T> &v,function<bool(T,T)>cmp=[](T x,T y){return x>y;}){
    int n=v.size();
    V<pii>ret(n,pii(0,n-1));
    stack<int>st;
    For(i,n){
        while(st.size()&&cmp(v[i],v[st.top()]))ret[st.top()].se=i-1,st.pop();
        if(st.size())ret[i].fi=st.top()+1;
        st.push(i);
    }
    return ret;
}
template<class T>
inline V<pii> cart_son(const V<T> &v,function<bool(T,T)>cmp=[](T x,T y){return x>y;}){
    int n=v.size();
    V<pii>ret(n,pii(-1,n));
    stack<int>st;
    For(i,n){
        while(st.size()&&cmp(v[i],v[st.top()]))ret[i].fi=st.top(),st.pop();
        if(st.size())ret[st.top()].se=i;
        st.push(i);
    }
    return ret;
}

struct lca_table{
    int n;
    V<int>dep,fa,siz,son,top;
    V<V<int>>to;
    inline lca_table(int n=0):n(n),dep(n),fa(n),siz(n),son(n,-1),top(n,-1),to(n){}
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
        to[x].pb(y),to[y].pb(x);
    }
    inline lca_table(const V<V<int>>&to,int rt=0):n(to.size()),dep(n),fa(n),siz(n),son(n,-1),top(n,-1),to(to){init(rt);}
    inline void init(int rt=0){
        function<void(int,int)>dfs1=[&](int p,int f){
            if(~f)dep[p]=dep[f]+1;
            fa[p]=f,siz[p]=1;
            for(int i:to[p])
                if(i!=f){
                    dfs1(i,p);
                    siz[p]+=siz[i];
                    if(!~son[p]||siz[i]>siz[son[p]])son[p]=i;
                }
        };
        dfs1(rt,-1);
        function<void(int,int)>dfs2=[&](int p,int k){
            top[p]=k;
            if(~son[p]){
                dfs2(son[p],k);
                for(int i:to[p])
                    if(!~top[i])
                        dfs2(i,i);
            }
        };
        dfs2(rt,rt);
    }
    inline int lca(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n);
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])swap(x,y);
            x=fa[top[x]];
        }
        return dep[x]<dep[y]?x:y;
    }
};

template<template<class,class>class T=RMQ>
struct lca_table_o1{
    int n;
    V<int>dfn;
    T<int,function<bool(int,int)>>rmq;
    V<V<int>>to;
    inline lca_table_o1(int n=0):n(n),dfn(n),to(n){}
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
        to[x].pb(y),to[y].pb(x);
    }
    inline lca_table_o1(const V<V<int>> &to,int rt=0):n(to.size()),dfn(n),to(to){init(rt);}
    inline void init(int rt=0){
        V<int>a;
        a.reserve(n-1);
        V<int>(n).swap(dfn);
        int cnt=0;
        auto dfs=[&](auto &&self,int p,int fa)->void{
            if(p!=rt)a.pb(fa);
            dfn[p]=cnt++;
            for(int i:to[p])if(i!=fa)self(self,i,p);
        };
        dfs(dfs,rt,-1);
        rmq=decltype(rmq)(a,[&](int x,int y){return dfn[x]<dfn[y];});
    }
    inline int lca(int x,int y)const{
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n);
        if(x==y)return x;
        x=dfn[x],y=dfn[y];
        if(x>y)swap(x,y);
        return rmq.query(x,y-1);
    }
};

struct tree_chain{
    int n;
    V<int>dep,dfn,fa,rev,siz,son,top;
    V<V<int>>to;
    inline tree_chain(int n=0):n(n),dep(n),dfn(n),fa(n),rev(n),siz(n),son(n,-1),top(n,-1),to(n){}
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
        to[x].pb(y),to[y].pb(x);
    }
    inline tree_chain(const V<V<int>>&to,int rt=0):n(to.size()),dep(n),dfn(n),fa(n),rev(n),siz(n),son(n,-1),top(n,-1),to(to){init(rt);}
    inline void init(int rt=0){
        function<void(int,int)>dfs1=[&](int p,int f){
            if(~f)dep[p]=dep[f]+1;
            fa[p]=f,siz[p]=1;
            for(int i:to[p])
                if(i!=f){
                    dfs1(i,p);
                    siz[p]+=siz[i];
                    if(!~son[p]||siz[i]>siz[son[p]])son[p]=i;
                }
        };
        dfs1(rt,-1);
        int cnt=0;
        function<void(int,int)>dfs2=[&](int p,int k){
            dfn[p]=cnt,rev[cnt++]=p,top[p]=k;
            if(~son[p]){
                dfs2(son[p],k);
                for(int i:to[p])
                    if(!~top[i])
                        dfs2(i,i);
            }
        };
        dfs2(rt,rt);
    }
    inline int lca(int x,int y)const{
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n);
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])swap(x,y);
            x=fa[top[x]];
        }
        return dep[x]<dep[y]?x:y;
    }
    inline int kthac(int p,int k){
        assert(0<=p),assert(p<n),assert(k>=0),assert(k<=dep[p]);
        while(k>dep[p]-dep[top[p]]){
            k-=dep[p]-dep[top[p]]+1;
            p=fa[top[p]];
        }
        return rev[dfn[p]-k];
    }
    inline V<pii> path(int x,int y,bool dir=0,bool lca=1){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n);
        V<pii>ret,ter;
        bool rv=0;
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])rv^=1,swap(x,y);
            if(dir){
                if(rv)ter.eb(dfn[top[x]],dfn[x]);
                else ret.eb(dfn[x],dfn[top[x]]);
            }
            else (rv?ter:ret).eb(dfn[top[x]],dfn[x]);
            x=fa[top[x]];
        }
        if(lca||x!=y){
            if(dep[x]>dep[y])rv^=1,swap(x,y);
            if(dir){
                if(rv)ter.eb(dfn[y],dfn[x]+!lca);
                else ret.eb(dfn[x]+!lca,dfn[y]);
            }
            else (rv?ret:ter).eb(dfn[x]+!lca,dfn[y]);
        }
        reverse(ALL(ter));
        ret.insert(ret.end(),ALL(ter));
        return ret;
    }
};

template<class T>
inline void virt_tree(V<int> &p,const T &t,V<V<int>> &to){
    sort(ALL(p),[&](int x,int y){return t.dfn[x]<t.dfn[y];});
    p.erase(unique(ALL(p)),p.end());
    auto add_edge=[&](int x,int y){to[x].pb(y),to[y].pb(x);};
    V<int>st;
    for(int i:p){
        if(st.size()){
            int anc=t.lca(i,st.back());
            if(anc!=st.back()){
                while(st.size()>1&&t.dfn[anc]<t.dfn[st[st.size()-2]])add_edge(st[st.size()-2],st.back()),st.qb();
                if(st.size()==1||t.dfn[anc]>t.dfn[st[st.size()-2]])V<int>().swap(to[anc]),add_edge(anc,st.back()),st.back()=anc;
                else add_edge(anc,st.back()),st.qb();
            }
        }
        V<int>().swap(to[i]),st.pb(i);
    }
    while(st.size()>1)add_edge(st[st.size()-2],st.back()),st.qb();
}

struct lca_table_forest{
    int n;
    V<int>dep,fa,id,siz,son,top;
    V<V<int>>to;
    inline lca_table_forest(int n=0):n(n),dep(n),fa(n,-1),id(n,-1),siz(n),son(n,-1),top(n,-1),to(n){}
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
        to[x].pb(y),to[y].pb(x);
    }
    inline lca_table_forest(const V<V<int>>&to):n(to.size()),dep(n),fa(n,-1),id(n,-1),siz(n),son(n,-1),top(n,-1),to(to){init();}
    inline void init(){ // 暂时不允许指定各子树的根
        function<void(int,int)>dfs1=[&](int p,int f){
            if(~f)dep[p]=dep[f]+1;
            fa[p]=f,siz[p]=1;
            for(int i:to[p])
                if(i!=f){
                    dfs1(i,p);
                    siz[p]+=siz[i];
                    if(!~son[p]||siz[i]>siz[son[p]])son[p]=i;
                }
        };
        For(i,n)if(!~fa[i])dfs1(i,-1);
        function<void(int,int)>dfs2=[&](int p,int k){
            id[p]=id[k],top[p]=k;
            if(~son[p]){
                dfs2(son[p],k);
                for(int i:to[p])
                    if(!~top[i])
                        id[i]=id[k],dfs2(i,i);
            }
        };
        For(i,n)if(!~top[i])id[i]=i,dfs2(i,i);
    }
    inline int lca(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(id[x]==id[y]);
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])swap(x,y);
            x=fa[top[x]];
        }
        return dep[x]<dep[y]?x:y;
    }
};

struct tree_chain_forest{
    int n;
    V<int>dep,dfn,fa,id,rev,siz,son,top;
    V<V<int>>to;
    inline tree_chain_forest(int n=0):n(n),dep(n),dfn(n),fa(n,-1),id(n,-1),rev(n),siz(n),son(n,-1),top(n,-1),to(n){}
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
        to[x].pb(y),to[y].pb(x);
    }
    inline tree_chain_forest(const V<V<int>>&to):n(to.size()),dep(n),dfn(n),fa(n,-1),id(n,-1),rev(n),siz(n),son(n,-1),top(n,-1),to(to){init();}
    inline void init(){ // 暂时不允许指定各子树的根
        function<void(int,int)>dfs1=[&](int p,int f){
            if(~f)dep[p]=dep[f]+1;
            fa[p]=f,siz[p]=1;
            for(int i:to[p])
                if(i!=f){
                    dfs1(i,p);
                    siz[p]+=siz[i];
                    if(!~son[p]||siz[i]>siz[son[p]])son[p]=i;
                }
        };
        For(i,n)if(!~fa[i])dfs1(i,-1);
        int cnt=0;
        function<void(int,int)>dfs2=[&](int p,int k){
            dfn[p]=cnt,rev[cnt++]=p;
            id[p]=id[k],top[p]=k;
            if(~son[p]){
                dfs2(son[p],k);
                for(int i:to[p])
                    if(!~top[i])
                        id[i]=id[k],dfs2(i,i);
            }
        };
        For(i,n)if(!~top[i])id[i]=i,dfs2(i,i);
    }
    inline int lca(int x,int y)const{
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(id[x]==id[y]);
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])swap(x,y);
            x=fa[top[x]];
        }
        return dep[x]<dep[y]?x:y;
    }
    inline int kthac(int p,int k){
        assert(0<=p),assert(p<n),assert(k>=0),assert(k<=dep[p]);
        while(k>dep[p]-dep[top[p]]){
            k-=dep[p]-dep[top[p]]+1;
            p=fa[top[p]];
        }
        return rev[dfn[p]-k];
    }
    inline V<pii> path(int x,int y,bool dir=0,bool lca=1){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(id[x]==id[y]);
        V<pii>ret,ter;
        bool rv=0;
        while(top[x]!=top[y]){
            if(dep[top[x]]<dep[top[y]])rv^=1,swap(x,y);
            if(dir){
                if(rv)ter.eb(dfn[top[x]],dfn[x]);
                else ret.eb(dfn[x],dfn[top[x]]);
            }
            else (rv?ter:ret).eb(dfn[top[x]],dfn[x]);
            x=fa[top[x]];
        }
        if(lca||x!=y){
            if(dep[x]>dep[y])rv^=1,swap(x,y);
            if(dir){
                if(rv)ter.eb(dfn[y],dfn[x]+!lca);
                else ret.eb(dfn[x]+!lca,dfn[y]);
            }
            else (rv?ret:ter).eb(dfn[x]+!lca,dfn[y]);
        }
        reverse(ALL(ter));
        ret.insert(ret.end(),ALL(ter));
        return ret;
    }
};

// root: n-1
inline V<int> pru2fa(const V<int> &_p){
    int n=_p.size()+2;
    V<int>deg(n),p=_p;p.pb(n-1);
    for(int i:_p)++deg[i];
    V<int>fa(n-1);
    int j=0;
    For(i,n-1){
        while(deg[j])++j;
        fa[j]=p[i];
        while(i<n-1&&!--deg[p[i]]&&p[i]<j){
            if(i+1<n-1)fa[p[i]]=p[i+1];
            ++i;
        }
        ++j;
    }
    return fa;
}
inline V<V<int>> pru2tr(const V<int> &p){
    int n=p.size()+2;
    V<int>fa=pru2fa(p);
    V<V<int>>to(n);
    For(i,n-1)to[i].pb(fa[i]),to[fa[i]].pb(i);
    return to;
}
inline V<int> fa2pru(const V<int> &fa){
    int n=fa.size()+1;
    V<int>deg(n);
    for(int i:fa)++deg[i];
    int j=0;
    V<int>p(n-2);
    For(i,n-2){
        while(deg[j])++j;
        p[i]=fa[j];
        while(i<n-2&&!--deg[p[i]]&&p[i]<j){
            if(i+1<n-2)p[i+1]=fa[p[i]];
            ++i;
        }
        ++j;
    }
    return p;
}
inline V<int> tr2pru(const V<V<int>> &to){
    int n=to.size();
    V<int>fa(n-1,-1);
    queue<int>q;
    q.push(n-1);
    while(q.size()){
        int p=q.front();q.pop();
        for(int i:to[p])if(i<n-1&&!~fa[i])fa[i]=p,q.push(i);
    }
    return fa2pru(fa);
}

inline V<int> ahu(int n,const V<V<int>> &to,int rt=0){
    assert(n>=0),assert(to.size()==n);
    For(i,n)for(int j:to[i])assert(0<=j),assert(j<n);
    V<int>a(n);
    static genID<V<int>,map<V<int>,int>>g;
    auto dfs=[&](auto &&self,int p,int fa)->void{
        V<int>tmp;
        for(int i:to[p])if(i!=fa){
            self(self,i,p);
            tmp.pb(a[i]);
        }
        sort(ALL(tmp));
        a[p]=g.get_id(tmp);
    };
    dfs(dfs,rt,-1);
    return a;
}
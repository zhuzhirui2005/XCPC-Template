template<class T,T e,T(*merge)(T,T)>
struct SGT{
    int n;
    V<T>t;
    inline SGT(int n=0):n(n),t(n<<2,e){}
    inline void push_up(int p){t[p]=merge(t[p<<1],t[p<<1|1]);}
    void build(int p,int l,int r,const V<T>&v){
        if(l==r){t[p]=v[l];return;}
        int mid=l+r>>1;
        build(p<<1,l,mid,v),build(p<<1|1,mid+1,r,v);
        push_up(p);
    }
    inline void build(const V<T>&v){
        assert(v.size()==n);
        build(1,0,n-1,v);
    }
    void build(int p,int l,int r){
        if(l==r){t[p]=e;return;}
        int mid=l+r>>1;
        build(p<<1,l,mid),build(p<<1|1,mid+1,r);
        push_up(p);
    }
    inline void build(){
        build(1,0,n-1);
    }
    void query(int p,int l,int r,int ql,int qr,T &ret){
        if(ql<=l&&r<=qr){ret=merge(ret,t[p]);return;}
        int mid=l+r>>1;
        if(ql<=mid)query(p<<1,l,mid,ql,qr,ret);
        if(qr>mid)query(p<<1|1,mid+1,r,ql,qr,ret);
    }
    inline T query(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<n);
        T ret=e;
        query(1,0,n-1,l,r,ret);
        return ret;
    }
    void modify(int p,int l,int r,int k,const T &v){
        if(l==r){t[p]=v;return;}
        int mid=l+r>>1;
        k<=mid?modify(p<<1,l,mid,k,v):modify(p<<1|1,mid+1,r,k,v);
        push_up(p);
    }
    inline void modify(int k,const T& v){
        assert(0<=k),assert(k<n);
        modify(1,0,n-1,k,v);
    }
};

template<class T,T e>
struct SGTlazy{
    int n;
    V<T>t,tag;
    inline SGTlazy(int n=0):n(n),t(n<<2,e),tag(n<<2,e){}
    inline void push_up(int p){}
    inline void add_tag(int p,const T &v){}
    inline void push_down(int p){add_tag(p<<1,tag[p]),add_tag(p<<1|1,tag[p]),tag[p]=e;}
    void build(int p,int l,int r,const V<T>&v){
        if(l==r){t[p]=v[l];return;}
        int mid=l+r>>1;
        build(p<<1,l,mid,v),build(p<<1|1,mid+1,r,v);
        push_up(p);
    }
    inline void build(const V<T>&v){
        assert(v.size()==n);
        build(1,0,n-1,v);
    }
    void query(int p,int l,int r,int ql,int qr,T &ret){
        if(ql<=l&&r<=qr){return;}
        push_down(p);
        int mid=l+r>>1;
        if(ql<=mid)query(p<<1,l,mid,ql,qr,ret);
        if(qr>mid)query(p<<1|1,mid+1,r,ql,qr,ret);
    }
    inline T query(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<n);
        T ret=e;
        query(1,0,n-1,l,r,ret);
        return ret;
    }
    void modify(int p,int l,int r,int k,const T &v){
        if(l==r){t[p]=v;return;}
        push_down(p);
        int mid=l+r>>1;
        k<=mid?modify(p<<1,l,mid,k,v):modify(p<<1|1,mid+1,r,k,v);
        push_up(p);
    }
    inline void modify(int k,const T& v){
        assert(0<=k),assert(k<n);
        modify(1,0,n-1,k,v);
    }
    void add(int p,int l,int r,int ql,int qr,const T &v){
        if(ql<=l&&r<=qr){add_tag(p,v);return;}
        push_down(p);
        int mid=l+r>>1;
        if(ql<=mid)add(p<<1,l,mid,ql,qr,v);
        if(qr>mid)add(p<<1|1,mid+1,r,ql,qr,v);
        push_up(p);
    }
    inline void add(int l,int r,const T &v){
        assert(0<=l),assert(l<=r),assert(r<n);
        add(1,0,n-1,l,r,v);
    }
    // int find_l(int p,int l,int r,int ql,int qr,const T &v){
    //     if(l==r)return l;
    //     push_down(p);
    //     int mid=l+r>>1;
    //     if(ql>mid||)return find_l(p<<1|1,mid+1,r,ql,qr,v);
    //     return find_l(p<<1,l,mid,ql,qr,v);
    // }
    // inline int find_l(int l,int r,const T &v){
    //     assert(0<=l),assert(l<=r),assert(r<n),assert(t[1]>=v);
    //     return find_l(1,0,n-1,l,r,v);
    // }
    // int find_r(int p,int l,int r,int ql,int qr,const T &v){
    //     if(l==r)return r;
    //     push_down(p);
    //     int mid=l+r>>1;
    //     if(qr<=mid||)return find_r(p<<1,l,mid,ql,qr,v);
    //     return find_r(p<<1|1,mid+1,r,ql,qr,v);
    // }
    // inline int find_r(int l,int r,const T &v){
    //     assert(0<=l),assert(l<=r),assert(r<n),assert(t[1]>=v);
    //     return find_r(1,0,n-1,l,r,v);
    // }
};

template<class T>
struct SGT_2n{
    int n;
    V<T>t,tag;
    inline int idx(int l,int r){return l+r|l!=r;}
    #define p idx(l,r)
    #define ls idx(l,mid)
    #define rs idx(mid+1,r)
    inline SGT_2n(int n=0):n(n),t(n<<1),tag(n<<1){}
    void build(int l,int r,const V<T>&v){
        if(l==r){t[p]=v[l];return;}
        int mid=l+r>>1;
        build(l,mid,v),build(mid+1,r,v);
        t[p]=max(t[ls],t[rs]);
    }
    inline void build(const V<T>&v){build(0,n-1,v);}
    void modify(int l,int r,int k,const T &v){
        if(l==r){t[p]=v;return;}
        int mid=l+r>>1;
        if(tag[p])t[ls]+=tag[p],t[rs]+=tag[p],tag[ls]+=tag[p],tag[rs]+=tag[p],tag[p]=0;
        k<=mid?modify(l,mid,k,v):modify(mid+1,r,k,v);
        t[p]=max(t[ls],t[rs]);
    }
    inline void modify(int k,const T& v){modify(0,n-1,k,v);}
    #undef p
    #undef ls
    #undef rs
};

template<class T>
struct ODT{
    map<T,T>p;
    inline ODT(){}
    inline ODT(const V<pair<T,T>> &v){
        for(const auto &[l,r]:v)insert(l,r);
    }
    inline V<pair<T,T>> insert(T l,T r){
        assert(l<=r);
        V<pair<T,T>>ret;
        T tmp=l-1;
        auto it=p.lower_bound(l-1);
        if(it!=p.end())ckmin(l,it->se);
        while(it!=p.end()&&it->se<=r+1){
            if(tmp+1<it->se)ret.eb(tmp+1,it->se-1);
            ckmax(r,tmp=it->fi);
            it=p.erase(it);
        }
        if(tmp<r)ret.eb(tmp+1,r);
        p.insert(it,{r,l});
        return ret;
    }
    inline V<pair<T,T>> erase(const T &l,const T &r){
        assert(l<=r);
        V<pair<T,T>>ret;
        auto it=p.lower_bound(l);
        if(it!=p.end()&&it->se<l)p.insert(it,{l-1,it->se});
        while(it!=p.end()&&it->fi<=r){
            ret.eb(max(it->se,l),it->fi);
            it=p.erase(it);
        }
        if(it!=p.end()&&it->se<=r){
            ret.eb(max(it->se,l),r);
            it->se=r+1;
        }
        return ret;
    }
    inline bool cover(const T &l,const T &r){
        assert(l<=r);
        auto it=p.lower_bound(r);
        return it!=p.end()&&it->se<=l;
    }
    inline bool intersect(const T &l,const T &r){
        assert(l<=r);
        auto it=p.lower_bound(l);
        return it!=p.end()&&it->se<=r;
    }
};

template<class T,T e,class U=less<T>>
struct lichao{
    U cmp;
    T l,r;
    struct LCT{
        LCT *ls,*rs;
        pair<T,T>v;
        inline LCT(T k=0,T b=e):ls(nullptr),rs(nullptr),v(k,b){}
    }*rt;
    inline lichao(T l,T r,const U &cmp=U()):cmp(cmp),l(l),r(r),rt(nullptr){}
    void insert(LCT *&p,T l,T r,pair<T,T> v){
        if(!p){
            p=new LCT(v.fi,v.se);
            return;
        }
        T mid=l+(r-l>>1);
        if(cmp(v.fi*mid+v.se,p->v.fi*mid+p->v.se))swap(v,p->v);
        if(l==r)return;
        if(cmp(v.fi*l+v.se,p->v.fi*l+p->v.se))insert(p->ls,l,mid,v);
        else if(cmp(v.fi*r+v.se,p->v.fi*r+p->v.se))insert(p->rs,mid+1,r,v);
    }
    inline void insert(const pair<T,T> &v){insert(rt,l,r,v);}
    T query(LCT *p,T l,T r,T x){
        if(!p)return e;
        T ret=p->v.fi*x+p->v.se;
        if(l==r)return ret;
        T mid=l+(r-l>>1);
        return min(ret,x<=mid?query(p->ls,l,mid,x):query(p->rs,mid+1,r,x),cmp);
    }
    inline T query(T x){return query(rt,l,r,x);}
};

template<class T>
struct dynamic_wavelet{
    T l,r;
    dynamic_wavelet *ls,*rs;
    V<int>pre;
    inline dynamic_wavelet(T l,T r):l(l),r(r),ls(nullptr),rs(nullptr){}
    dynamic_wavelet(T l,T r,typename V<T>::iterator ql,typename V<T>::iterator qr):dynamic_wavelet(l,r){
        if(l==r)return;
        T mid=l+r>>1;
        pre.reserve(qr-ql);
        for(auto it=ql;it!=qr;++it)pre.pb((pre.size()?pre.back():0)+(*it<=mid));
        auto qm=stable_partition(ql,qr,[&](T k){return k<=mid;}); // this will modify arg arr
        if(ql<qm)ls=new dynamic_wavelet(l,mid,ql,qm);
        if(qm<qr)rs=new dynamic_wavelet(mid+1,r,qm,qr);
    }
    void append(T k){
        if(l==r)return;
        T mid=l+r>>1;
        if(k<=mid){
            pre.pb(pre.size()?pre.back()+1:1);
            if(!ls)ls=new dynamic_wavelet(l,mid);
            ls->append(k);
        }
        else{
            pre.pb(pre.size()?pre.back():0);
            if(!rs)rs=new dynamic_wavelet(mid+1,r);
            rs->append(k);
        }
    }
    void pop(){
        if(pre.empty())return;
        int lst=pre.back();pre.qb();
        if((pre.size()?pre.back():0)<lst){
            ls->pop();
            if(ls->l<ls->r&&ls->pre.empty())ls=nullptr;
        }
        else{
            rs->pop();
            if(rs->l<rs->r&&rs->pre.empty())rs=nullptr;
        }
    }
    T kth(int ql,int qr,int k){ // 0-indexed
        if(l==r)return l;
        int cl=ql?pre[ql-1]:0,cr=pre[qr];
        return cr-cl>k?ls->kth(cl,cr-1,k):rs->kth(ql-cl,qr-cr,k-(cr-cl));
    }
    int count(int ql,int qr,T k){
        if(ql>qr||k<l)return 0;
        if(r<=k)return qr-ql+1;
        int cl=ql?pre[ql-1]:0,cr=pre[qr];
        return (ls?ls->count(cl,cr-1,k):0)+(rs?rs->count(ql-cl,qr-cr,k):0);
    }
    int count(int ql,int qr,T x,T y){return count(ql,qr,y)-count(ql,qr,x-1);}
};

template<class T,int n>
struct wavelet{
    static_assert(n>0);
    static_assert(n<=(is_same_v<T,int>?31:is_same_v<T,unsigned>?32:is_same_v<T,ll>?63:is_same_v<T,ull>?64:-1));
    array<int,n>p;
    array<V<pair<ull,int>>,n>t;
    inline wavelet(V<T> &v){
        for(T i:v)assert(i>=0);
        int m=v.size();
        Rep(i,n){
            t[i].resize((m>>6)+1);
            For(j,m)if(!(v[j]>>i&1))t[i][j>>6].fi|=1ull<<(j&63);
            FOR(j,1,(m>>6)+1)t[i][j].se=t[i][j-1].se+__builtin_popcountll(t[i][j-1].fi);
            p[i]=stable_partition(ALL(v),[&](T k){return !(k>>i&1);})-v.begin(); // this will modify arg arr
        }
    }
    inline int get(int i,int k){ // [0, k)
        return t[i][k>>6].se+__builtin_popcountll(t[i][k>>6].fi&((1ull<<(k&63))-1));
    } 
    inline T kth(int ql,int qr,int k){ // 0-indexed
        T ret=0;
        Rep(i,n){
            int cl=get(i,ql),cr=get(i,qr+1);
            if(cr-cl>k)ql=cl,qr=cr-1;
            else k-=cr-cl,ql+=p[i]-cl,qr+=p[i]-cr,ret|=T(1)<<i;
        }
        return ret;
    }
    inline int count(int ql,int qr,T k){
        int ret=0;
        Rep(i,n){
            if(ql>qr)break;
            int cl=get(i,ql),cr=get(i,qr+1);
            if(k>>i&1)ql+=p[i]-cl,qr+=p[i]-cr,ret+=cr-cl;
            else ql=cl,qr=cr-1;
        }
        return ret+(ql>qr?0:qr-ql+1);
    }
    inline int count(int ql,int qr,T x,T y){return count(ql,qr,y)-(x?count(ql,qr,x-1):0);}
};
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
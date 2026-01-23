template<class T>
inline void radix_sort(V<T> &v){
    int n=v.size();
    if(is_same_v<T,int>){
        V<unsigned>tmp(n);
        memcpy(tmp.data(),v.data(),n<<2);
        for(unsigned &i:tmp)i^=1u<<31;
        radix_sort(tmp);
        for(unsigned &i:tmp)i^=1u<<31;
        memcpy(v.data(),tmp.data(),n<<2);
    }
    else if(is_same_v<T,ll>){
        V<ull>tmp(n);
        memcpy(tmp.data(),v.data(),n<<3);
        for(ull &i:tmp)i^=1ull<<63;
        radix_sort(tmp);
        for(ull &i:tmp)i^=1ull<<63;
        memcpy(v.data(),tmp.data(),n<<3);
    }
    else{
        assert((is_same_v<T,unsigned>)||(is_same_v<T,ull>));
        static int cnt[1<<16];
        V<T>tmp(n);
        For(i,(is_same_v<T,ull>)?4:2){
            memset(cnt,0,1<<18);
            for(T j:v)++cnt[(j>>(i<<4))&65535];
            int pre=0;
            for(int &j:cnt)swap(j,pre),pre+=j;
            for(T j:v)tmp[cnt[(j>>(i<<4))&65535]++]=j;
            v.swap(tmp);
        }
    }
}

template<class T>
inline V<V<T>> rot(const V<V<T>> &v){
    V<V<T>>ret(v[0].size(),V<T>(v.size()));
    For(i,v.size())
        For(j,v[0].size())
            ret[j][v.size()-i-1]=v[i][j];
    return ret;
}
template<class T>
inline V<V<T>> rot(const V<V<T>> &v,int k){
    assert(k>=0);
    k&=3;
    auto ret=v;
    while(k--)ret=rot(ret);
    return ret;
}

struct ListNode {
     int val;
     ListNode *next;
     ListNode() : val(0), next(nullptr) {}
     ListNode(int x) : val(x), next(nullptr) {}
     ListNode(int x, ListNode *next) : val(x), next(next) {}
};
inline V<int> ltov(ListNode *hd){
    V<int>ret;
    while(hd)ret.pb(hd->val),hd=hd->next;
    return ret;
}
inline ListNode* vtol(const V<int>& v){
    if(v.empty())return NULL;
    ListNode *hd=new ListNode(v[0]),*p=hd;
    FOR(i,1,v.size())p->next=new ListNode(v[i]),p=p->next;
    return hd;
}

template<class T>
inline T cantor(const V<int> &v){
    int d=*min_element(ALL(v)),n=v.size();
    V<bool>vis(n);
    for(int i:v)vis[i-d]=true;
    For(i,n)if(!vis[i])return -1;
    V<T>fac(n);
    fac[0]=1;
    BIT3<int>t(n);
    FOR(i,1,n+1){
        if(i<n)fac[i]=fac[i-1]*i;
        ++t.c[i];
        if(i+(i&-i)<=n)t.c[i+(i&-i)]+=t.c[i];
    }
    T ret=0;
    For(i,n){
        t.add(v[i]-d,-1);
        ret+=fac[n-i-1]*t.query(v[i]-d);
    }
    return ret;
}

inline V<int> decantor(int n,ll k){
    V<ll>fac(n+1);
    fac[0]=1;
    FOR(i,1,n+1)fac[i]=fac[i-1]*i;
    if(k>=fac[n])return {-1};
    V<int>ret(n);
    V<bool>vis(n);
    For(i,n){
        int d=k/fac[n-i-1]+1,j=-1;
        k%=fac[n-i-1];
        do d-=!vis[++j];while(d);
        ret[i]=j,vis[j]=true;
    }
    return ret;
}
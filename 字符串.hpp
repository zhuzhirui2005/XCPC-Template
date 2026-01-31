struct manacher{
    int n;
    V<int>p;
    inline manacher(const string &s):n(s.size()),p(n<<1|1,1){
        string t(n<<1|1,'#');
        For(i,n)t[i<<1|1]=s[i];
        for(int i=0,mid=-1,mx=-1;i<p.size();++i){
            if(i<=mx)p[i]=min(p[(mid<<1)-i],mx-i)+1;
            while(i>=p[i]&&i+p[i]<p.size()&&t[i-p[i]]==t[i+p[i]])++p[i];
            if(i+--p[i]>mx)mid=i,mx=i+p[i];
        }
    }
    inline int odd(int k){
        assert(0<=k&&k<n);
        return p[k<<1|1];
    }
    inline int even(int k){
        assert(0<=k&&k+1<n);
        return p[k+1<<1];
    }
    inline bool isp(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<n);
        return p[l+r+1]>=r-l+1;
    }
};

inline V<int> get_kmp(const string &s){
    int n=s.size();
    V<int>kmp(n);
    for(int i=1,j=0;i<n;++i){
        while(j&&s[j]!=s[i])j=kmp[j-1];
        if(s[j]==s[i])++j;
        kmp[i]=j;
    }
    return kmp;
}
inline V<int> find_kmp(const V<int> &kmp,const string &s,const string &t){
    int n=s.size(),m=t.size();
    V<int>ret;
    for(int i=0,j=0;i<n;++i){
        while(j&&t[j]!=s[i])j=kmp[j-1];
        if(t[j]==s[i])++j;
        if(j==m)ret.pb(i);
    }
    return ret;
}

inline V<int> get_z(const string &s){
    int n=s.size();
    V<int>z(n);
    z[0]=n;
    for(int i=1,l=0,r=-1;i<n;++i){
        if(i<=r)z[i]=min(r-i+1,z[i-l]);
        while(i+z[i]<n&&s[z[i]]==s[i+z[i]])++z[i];
        if(ckmax(r,i+z[i]-1))l=i;
    }
    return z;
}
inline V<int> get_ext(const V<int> &z,const string &s,const string &t){
    int n=s.size(),m=t.size();
    V<int>ext(n);
    for(int i=0,l=0,r=-1;i<n;++i){
        if(i<=r)ext[i]=min(r-i+1,z[i-l]);
        while(i+ext[i]<n&&ext[i]<m&&t[ext[i]]==s[i+ext[i]])++ext[i];
        if(ckmax(r,i+ext[i]-1))l=i;
    }
    return ext;
}

inline pair<V<int>,V<int>> suffix_array(const string &s,int m=128){
    int n=s.size();
    V<int>cnt(m),id(n),rk(n),sa(n),tmp(n);
    if(n==1)return {rk,sa};
    For(i,n)++cnt[rk[i]=s[i]];
    FOR(i,1,m)cnt[i]+=cnt[i-1];
    Rep(i,n)sa[--cnt[rk[i]]]=i;
    for(int p=0,w=1;p+1<n;m=p+1,w<<=1){
        id.clear();
        FOR(i,n-w,n)id.pb(i);
        For(i,n)if(sa[i]>=w)id.pb(sa[i]-w);
        cnt.assign(m,0);
        For(i,n)++cnt[rk[i]];
        FOR(i,1,m)cnt[i]+=cnt[i-1];
        Rep(i,n)sa[--cnt[rk[id[i]]]]=id[i];
        rk.swap(tmp);
        p=rk[sa[0]]=0;
        FOR(i,1,n){
            if(tmp[sa[i]]==tmp[sa[i-1]]&&(sa[i]+w<n?tmp[sa[i]+w]:-1)==(sa[i-1]+w<n?tmp[sa[i-1]+w]:-1))rk[sa[i]]=p;
            else rk[sa[i]]=++p;
        }
    }
    return {rk,sa};
}
inline V<int> get_ht(const string &s,const V<int> &rk,const V<int> &sa){
    int n=s.size();
    V<int>ht(n-1);
    for(int i=0,k=0;i<n;++i){
        if(k)--k;
        int r=rk[i];
        if(!r)continue;
        int j=sa[r-1];
        while(max(i,j)+k<n&&s[i+k]==s[j+k])++k;
        ht[r-1]=k;
    }
    return ht;
}

template<int l=1>
inline V<int> shift_and(const string s,const string &t,char wc='*'){
    for(char c:s)assert(c==wc||islower(c));
    for(char c:t)assert(c==wc||islower(c));
    int n=s.size(),m=t.size();
    static const int lim=1'000'000;
    if(l<=lim&&l<n)return shift_and<l<<1>(s,t);
    bitset<l>a;
    V<bitset<l>>b(26);
    For(i,n){
        if(s[i]==wc)For(j,26)b[j].set(i);
        else b[s[i]-'a'].set(i);
    }
    a.flip();
    For(i,m)if(t[i]!=wc)a&=b[t[i]-'a']>>i;
    V<int>ret;
    For(i,n-m+1)if(a[i])ret.pb(i);
    return ret;
}

inline V<int> ntt_match(const string &s,const string &t,int g,char wc='*'){
    int n=s.size(),m=t.size();
    V<mi>a(n-m+1),b(n),c(m),d;
    For(i,n){
        if(s[i]==wc)b[i]=0;
        else b[i]=(mi)s[i]*s[i]*s[i];
    }
    For(i,m){
        if(t[i]==wc)c[i]=0;
        else c[i]=t[i];
    }
    d=poly_conv_sub(b,c,g);
    For(i,n-m+1)a[i]+=d[i];
    For(i,n){
        if(s[i]==wc)b[i]=0;
        else b[i]=(mi)s[i]*s[i];
    }
    For(i,m){
        if(t[i]==wc)c[i]=0;
        else c[i]=(mi)t[i]*t[i];
    }
    d=poly_conv_sub(b,c,g);
    For(i,n-m+1)a[i]-=d[i]+d[i];
    For(i,n){
        if(s[i]==wc)b[i]=0;
        else b[i]=s[i];
    }
    For(i,m){
        if(t[i]==wc)c[i]=0;
        else c[i]=(mi)t[i]*t[i]*t[i];
    }
    d=poly_conv_sub(b,c,g);
    For(i,n-m+1)a[i]+=d[i];
    V<int>ret;
    For(i,n-m+1)if(!a[i])ret.pb(i);
    return ret;
}

struct eertree{
    V<int>d,fail,len,link;
    int lst;
    V<array<int,26>>nxt;
    string s;
    inline eertree():d(2),fail(2),len{-1,0},link(2),lst(1),nxt{{},{}}{}
    inline bool extend(char c){
        int pos=s.size();
        s+=c,c-='a';
        while(pos<len[lst]+1||s[pos-1-len[lst]]!=s[pos])lst=fail[lst];
        if(!nxt[lst][c]){
            len.pb(len[lst]+2),nxt[lst][c]=nxt.size();
            if(len.back()>1){
                for(lst=fail[lst];pos<len[lst]+1||s[pos-1-len[lst]]!=s[pos];lst=fail[lst]);
                fail.pb(nxt[lst][c]);
            }
            else fail.pb(1);
            d.pb(len.back()-len[fail.back()]),link.pb(d.back()>d[fail.back()]?fail.back():link[fail.back()]),lst=nxt.size(),nxt.pb({});
            return true;
        }
        else{
            lst=nxt[lst][c];
            return false;
        }
    }
    inline eertree(const string &t):eertree(){for(char c:t)extend(c);}
};
template<class T>
inline V<int> minpal(const T &s){
    eertree e;
    V<int>f{0},g(2);
    int n=s.size();
    f.reserve(n);
    for(char c:s){
        int i=f.size();f.pb(inf);
        if(e.extend(c))g.pb(inf);
        for(int j=e.lst;j>1;j=e.link[j]){
            g[j]=f[i-e.d[j]-e.len[e.link[j]]];
            if(e.fail[j]!=e.link[j])ckmin(g[j],g[e.fail[j]]);
            ckmin(f.back(),g[j]+1);
        }
    }
    return f;
}

template<template<class> class T=RMQ>
inline V<array<int,3>> runs(const string &s,int m=128){
    int n=s.size();
    string r=s,t=s;
    reverse(ALL(r));
    for(char &c:t)c=m-1-c;
    auto [rkr,sar]=suffix_array(r,m);
    auto [rks,sas]=suffix_array(s,m);
    auto [rkt,sat]=suffix_array(t,m);
    V<int>htr=get_ht(r,rkr,sar),hts=get_ht(s,rks,sas);
    T<int>rmqr(htr),rmqs(hts);
    auto lcp=[&](int x,int y){
        if(x==y)return n-x;
        x=rks[x],y=rks[y];
        if(x>y)swap(x,y);
        return rmqs.query(x,y-1);
    };
    auto lcs=[&](int x,int y){
        if(x==y)return x+1;
        x=rkr[n-1-x],y=rkr[n-1-y];
        if(x>y)swap(x,y);
        return rmqr.query(x,y-1);
    };
    V<array<int,3>>ret;
    auto calc=[&](const V<int> &rk){
        V<int>st;
        Rep(i,n){
            while(st.size()&&rk[i]<rk[st.back()])st.qb();
            if(st.size()){
                int j=st.back(),l=lcs(i,j),r=lcp(i,j);
                if(l+r>j-i)ret.pb({i-l+1,j+r-1,j-i});
            }
            st.pb(i);
        }
    };
    calc(rks),calc(rkt);
    sort(ALL(ret));
    ret.erase(unique(ALL(ret)),ret.end());
    return ret;
}

template<class T>
inline V<ull> lcsarr(T al,T ar,T bl,T br){
    if(al==ar||bl==br)return {};
    int n=ar-al,k=n+63>>6;
    V<ull>f(k);
    char mn=*min_element(al,ar),mx=*max_element(al,ar);
    V g(mx-mn+1,V<ull>(k));
    for(auto it=al;it!=ar;++it)g[*it-mn][it-al>>6]|=1ull<<((it-al)&63);
    for(auto it=bl;it!=br;++it)if(mn<=*it&&*it<=mx){
        char i=*it-mn;
        bool w=0,z=1;
        For(j,k){
            ull x=f[j]<<1|z,y=f[j]|g[i][j],v=y-x-w;
            z=f[j]>>63,f[j]=y&~v,w=(w&(x==ULLONG_MAX))|(y<v);
        }
    }
    return f;
}

inline int lcs(const string &a,const string &b){
    V<ull>f=lcsarr(ALL(a),ALL(b));
    return accumulate(ALL(f),0,[&](int x,ull y){return x+__builtin_popcountll(y);});
}

inline string hirschberg(string::iterator al,string::iterator ar,string::iterator bl,string::iterator br){
    int n=ar-al,m=br-bl;
    if(min(n,m)<1)return "";
    if(n==1){
        for(auto it=bl;it!=br;++it)if(*it==*al)return string(1,*it);
        return "";
    }
    if(m==1){
        for(auto it=al;it!=ar;++it)if(*it==*bl)return string(1,*it);
        return "";
    }
    string::iterator am=al+(n>>1);
    #define mri make_reverse_iterator
    V<ull>fl=lcsarr(bl,br,al,am),fr=lcsarr(mri(br),mri(bl),mri(ar),mri(am));
    #undef mri
    int cnt=accumulate(ALL(fr),0,[&](int x,ull y){return x+__builtin_popcountll(y);}),mx=cnt,pos=0;
    For(i,m){
        if(fl[i>>6]>>(i&63)&1)++cnt;
        if(fr[m-1-i>>6]>>((m-1-i)&63)&1)--cnt;
        if(ckmax(mx,cnt))pos=i+1;
    }
    string::iterator bm=bl+pos;
    return hirschberg(al,am,bl,bm)+hirschberg(am,ar,bm,br);
}

inline int lcs63(const string &a,const string &b){
    if(a.empty()||b.empty())return 0;
    int n=a.size(),k=(n+62)/63;
    V<ull>f(k);
    char mn=*min_element(ALL(a)),mx=*max_element(ALL(a));
    V g(mx-mn+1,V<ull>(k));
    For(i,n)g[a[i]-mn][i/63]|=1ull<<i%63;
    for(char i:b)if(mn<=i&&i<=mx){
        i-=mn;
        bool z=1;
        For(j,k){
            ull x=f[j],y=f[j]|g[i][j];
            ((x<<=1)|=z)+=(~y)&((1ull<<63)-1);
            f[j]=x&y,z=x>>63;
        }
    }
    return accumulate(ALL(f),0,[&](int x,ull y){return x+__builtin_popcountll(y);});
}

struct subseq_table{
    V<V<int>>nxt;
    inline subseq_table(const string &v):nxt(128){
        For(i,v.size())nxt[v[i]].pb(i);
    }
    inline int lcp(const string &v){
        int nw=0,ret=0;
        for(char i:v){
            assert(i>=0&&i<128);
            auto it=lower_bound(ALL(nxt[i]),nw);
            if(it==nxt[i].end())break;
            nw=*it+1,++ret;
        }
        return ret;
    }
    inline bool query(const string &v){
        return lcp(v)==v.size();
    }
};

template<class T,class container=unordered_map<T,int>>
struct subseq_Table{
    genID<T,container>g;
    V<V<int>>nxt;
    inline subseq_Table(const V<T> &v){
        For(i,v.size()){
            int k=g.get_id(v[i]);
            if(k==nxt.size())nxt.pb({});
            nxt[k].pb(i);
        }
    }
    inline int lcp(const V<T> &v){
        int nw=0,ret=0;
        for(const T &i:v){
            if(!g.count(i))break;
            int k=g.get_id(i);
            auto it=lower_bound(ALL(nxt[k]),nw);
            if(it==nxt[k].end())break;
            nw=*it+1,++ret;
        }
        return ret;
    }
    inline bool query(const V<T> &v){
        return lcp(v)==v.size();
    }
};
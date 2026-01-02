/*
p = r * 2^k + 1£¬g ÎªÔ­¸ù¡£

p                   r   k   g
3                   1   1   2
5                   1   2   2
17                  1   4   3
97                  3   5   5
193                 3   6   5
257                 1   8   3
7681                15  9   17
12289               3   12  11
40961               5   13  3
65537               1   16  3
786433              3   18  10
5767169             11  19  3
7340033             7   20  3
23068673            11  21  3
104857601           25  22  3
167772161           5   25  3
469762049           7   26  3
998244353           119 23  3
1004535809          479 21  3
2013265921          15  27  31
2281701377          17  27  3
3221225473          3   30  5
75161927681         35  31  3
77309411329         9   33  7
206158430209        3   36  22
2061584302081       15  37  7
2748779069441       5   39  3
6597069766657       3   41  5
39582418599937      9   42  5
79164837199873      9   43  5
263882790666241     15  44  7
1231453023109121    35  45  3
1337006139375617    19  46  3
3799912185593857    27  47  5
4222124650659841    15  48  19
7881299347898369    7   50  6
31525197391593473   7   52  3
180143985094819841  5   55  6
1945555039024054273 27  56  5
4179340454199820289 29  57  3
*/

inline V<mi> poly_conv_add(V<mi> a,V<mi> b,int g){ // c[k]=¡Æ(a[i]*b[j]) for i+j=k verified with lg3803
    int n=a.size(),m=b.size();
    assert(n&&m);
    if(max(n,m)<17){
        V<mi>c(n+m-1);
        For(i,n)For(j,m)c[i+j]+=a[i]*b[j];
        return c;
    }
    int lg=0,nm=1;
    while(nm<n+m-1)++lg,nm<<=1;
    a.resize(nm),b.resize(nm);
    static V<V<int>>btf;
    while(btf.size()<=lg){
        int l=1<<btf.size();
        btf.pb({});
        V<int>&bf=btf.back();
        bf.resize(l);
        For(i,l)bf[i]=(bf[i>>1]>>1)|((i&1)?l>>1:0);
    }
    const V<int> &bf=btf[lg];
    auto NTT=[&](V<mi> &f,mi coef){
        For(i,nm)if(i<bf[i])swap(f[i],f[bf[i]]);
        for(int k=1,l=2;k<nm;k<<=1,l<<=1){
            mi wn=coef^((mod-1)/l);
            for(int i=0;i<nm;i+=l){
                mi w=1;
                For(j,k){
                    mi x=f[i|j],y=w*f[i|j|k];
                    f[i|j]=x+y,f[i|j|k]=x-y;
                    w*=wn;
                }
            }
        }
    };
    NTT(a,g),NTT(b,g);
    For(i,nm)a[i]*=b[i];
    NTT(a,mi(1)/g);
    a.resize(n+m-1);
    mi invnm=mi(1)/nm;
    for(mi &i:a)i*=invnm;
    return a;
}

inline V<mi> poly_conv_sub(const V<mi> &a,V<mi> b,int g){ // c[k]=¡Æ(a[i]*b[j]) for i-j=k verified with gym105386H
    int n=a.size(),m=b.size();
    assert(n&&m);
    if(n<m)return {};
    reverse(ALL(b));
    b=poly_conv_add(a,b,g);
    // (-b.size(),a.size()) -> [0,a.size())
    b.erase(b.begin(),b.begin()+m-1);
    return b;
}

inline int find_g(int m){
    auto phi=[&](int k){
        int ret=k;
        for(int i=2;i*i<=k;++i)if(k%i==0){ret-=ret/i;do k/=i;while(k%i==0);}
        if(k>1)ret-=ret/k;
        return ret;
    };
    int p=phi(m);
    V<int>fac;
    {
        int j=p;
        for(int i=2;i*i<=j;++i)if(j%i==0){fac.pb(p/i);do j/=i;while(j%i==0);}
        if(j>1)fac.pb(p/j);
    }
    auto check_g=[&](int g){
        auto qpow=[&](int x,int y){
            int z=1;
            for(;y;x=1ll*x*x%m,y>>=1)if(y&1)z=1ll*z*x%m;
            return z;
        };
        if(qpow(g,p)!=1)return false;
        for(int i:fac)if(qpow(g,i)==1)return false;
        return true;
    };
    FOR(i,1,m)if(check_g(i))return i;
    return -1;
}
inline V<mi> poly_conv_mul(V<mi> a,V<mi> b,int g,int p,int pg=-1){ // c[k]=¡Æ(a[i]*b[j]) for i*j%p=k (p should be prime) verified by qoj9247
    int n=a.size(),m=b.size();
    assert(n&&m),assert(p>1);
    for(int i=2;i*i<=p;++i)assert(p%i);
    if(!~pg)pg=find_g(p);
    assert(~pg);
    V<int>exp(p-1),lg(p);
    lg[0]=-1;
    for(int i=1,j=0;j<p-1;i=1ll*i*pg%p,++j)exp[j]=i,lg[i]=j;
    if(n>p){
        FOR(i,p,n)a[i%p]+=a[i];
        a.resize(p),n=p;
    }
    if(m>p){
        FOR(i,p,m)b[i%p]+=b[i];
        b.resize(p),m=p;
    }
    V<mi>_a(p-1),_b(p-1);
    FOR(i,1,n)_a[lg[i]]=a[i];
    FOR(i,1,m)_b[lg[i]]=b[i];
    V<mi>c=poly_conv_add(_a,_b,g);
    FOR(i,p-1,c.size())c[i-(p-1)]+=c[i];
    V<mi>d(p);
    d[0]=a[0]*reduce(ALL(b))+b[0]*reduce(ALL(a))-a[0]*b[0];
    For(i,p-1)d[exp[i]]=c[i];
    return d;
}

inline V<mi> poly_conv_div(V<mi> a,V<mi> b,int g,int p,int pg=-1){ // c[k]=¡Æ(a[i]*b[j]) for i/j%p=k (p should be prime) not verified
    int n=a.size(),m=b.size();
    assert(n&&m),assert(p>1);
    for(int i=2;i*i<=p;++i)assert(p%i);
    if(n>p){
        FOR(i,p,n)a[i%p]+=a[i];
        a.resize(p),n=p;
    }
    if(m>p){
        FOR(i,p,m)b[i%p]+=b[i];
        b.resize(p),m=p;
    }
    assert(!b[0]);
    static V<int>inv{0,1};
    while(inv.size()<p)inv.pb(1ll*(p-p/inv.size())*inv[p%inv.size()]%p);
    V<mi>_b(p);
    FOR(i,1,m)_b[inv[i]]=b[i];
    return poly_conv_mul(a,_b,g,p,pg);
}

inline V<mi> poly_conv_and(V<mi> a,V<mi> b){ // c[k]=¡Æ(a[i]*b[j]) for i&j=k verified with lg4717
    int n=a.size(),m=b.size();
    assert(n&&m);
    int nm=1;
    while(nm<max(n,m))nm<<=1;
    a.resize(nm),b.resize(nm);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<nm;k<<=1,l<<=1)for(int i=0;i<nm;i+=l)For(j,k)f[i|j]+=f[i|j|k]*coef;
    };
    FWT(a,1),FWT(b,1);
    For(i,nm)a[i]*=b[i];
    FWT(a,mod-1);
    return a;
}

inline V<mi> poly_conv_or(V<mi> a,V<mi> b){ // c[k]=¡Æ(a[i]*b[j]) for i|j=k verified with lg4717
    int n=a.size(),m=b.size();
    assert(n&&m);
    int nm=1;
    while(nm<max(n,m))nm<<=1;
    a.resize(nm),b.resize(nm);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<nm;k<<=1,l<<=1)for(int i=0;i<nm;i+=l)For(j,k)f[i|j|k]+=f[i|j]*coef;
    };
    FWT(a,1),FWT(b,1);
    For(i,nm)a[i]*=b[i];
    FWT(a,mod-1);
    return a;
}

inline V<mi> poly_conv_xor(V<mi> a,V<mi> b){ // c[k]=¡Æ(a[i]*b[j]) for i^j=k verified with lg4717
    int n=a.size(),m=b.size();
    assert(n&&m);
    int nm=1;
    while(nm<max(n,m))nm<<=1;
    a.resize(nm),b.resize(nm);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<nm;k<<=1,l<<=1)for(int i=0;i<nm;i+=l)For(j,k){
            mi x=f[i|j],y=f[i|j|k];
            f[i|j]=(x+y)*coef,f[i|j|k]=(x-y)*coef;
        }
    };
    FWT(a,1),FWT(b,1);
    For(i,nm)a[i]*=b[i];
    FWT(a,mod+1>>1);
    return a;
}

inline V<mi> poly_conv_gcd(const V<mi> &a,const V<mi> &b){ // c[k]=¡Æ(a[i]*b[j]) for gcd(i,j)=k verified with lc418t4
    int n=a.size(),m=b.size();
    assert(n&&m);
    int nm=max(n,m);
    V<mi>_a=a,_b=b;
    _a.resize(nm),_b.resize(nm);
    V<int>pri;
    V<bool>vis(nm);
    FOR(i,2,nm)if(!vis[i]){
        pri.pb(i);
        for(int k=(nm-1)/i,j=k*i;k;j-=i,--k)_a[k]+=_a[j],_b[k]+=_b[j],vis[j]=true;
    }
    FOR(i,1,nm)_a[i]*=_b[i];
    for(int i:pri)for(int j=i,k=1;j<nm;j+=i,++k)_a[k]-=_a[j];
    _a[0]=a[0]*b[0];
    FOR(i,1,nm){
        if(i<n)_a[i]+=b[0]*a[i];
        if(i<m)_a[i]+=a[0]*b[i];
    }
    return _a;
}

inline V<mi> poly_conv_lcm(const V<mi> &a,const V<mi> &b){ // c[k]=¡Æ(a[i]*b[j]) for lcm(i,j)=k not verified
    int n=a.size(),m=b.size();
    assert(n&&m);
    int nm=max(n,m);
    V<mi>_a=a,_b=b;
    _a.resize(nm),_b.resize(nm);
    V<int>pri;
    V<bool>vis(nm);
    FOR(i,2,nm)if(!vis[i]){
        pri.pb(i);
        for(int j=i,k=1;j<nm;j+=i,++k)_a[j]+=_a[k],_b[j]+=_b[k],vis[j]=true;
    }
    FOR(i,1,nm)_a[i]*=_b[i];
    for(int i:pri)for(int k=(nm-1)/i,j=k*i;k;j-=i,--k)_a[j]-=_a[k];
    _a[0]=a[0]*reduce(ALL(b))+b[0]*reduce(ALL(a))-a[0]*b[0];
    return _a;
}

inline V<mi> poly_inv(const V<mi> &a,int g){ // b=1/a verified with lg4238
    assert(a.size()),assert(a[0].val);
    V<mi>b{1/a[0]};
    mi invg=mi(1)/g,invm=1;
    int m=1;
    while(b.size()<a.size()){
        int n=min(a.size(),b.size()<<1);
        while(m<=n-1<<1)invm*=mod+1>>1,m<<=1;
        V<mi>c(a.begin(),a.begin()+n);
        b.resize(m),c.resize(m);
        V<int>bf(m);
        For(i,m)bf[i]=(bf[i>>1]>>1)|((i&1)?m>>1:0);
        auto NTT=[&](V<mi> &f,mi coef){
            For(i,m)if(i<bf[i])swap(f[i],f[bf[i]]);
            for(int k=1,l=2;k<m;k<<=1,l<<=1){
                mi wn=coef^((mod-1)/l);
                for(int i=0;i<m;i+=l){
                    mi w=1;
                    For(j,k){
                        mi x=f[i|j],y=w*f[i|j|k];
                        f[i|j]=x+y,f[i|j|k]=x-y;
                        w*=wn;
                    }
                }
            }
        };
        NTT(b,g),NTT(c,g);
        For(i,m)b[i]*=2-b[i]*c[i];
        NTT(b,invg);
        b.resize(n);
        for(mi &i:b)i*=invm;
    }
    return b;
}

inline V<mi> poly_diff(const V<mi> &a){ // b=a'
    int n=a.size();
    assert(n);
    if(n==1)return {0};
    V<mi>b(n-1);
    For(i,n-1)b[i]=a[i+1]*(i+1);
    return b;
}

inline V<mi> poly_intg(const V<mi> &a){ // b=¡Òa
    int n=a.size();
    assert(n);
    static V<mi>inv{0,1};
    while(inv.size()<=n)inv.pb((mod-mod/inv.size())*inv[mod%inv.size()]);
    V<mi>b(n+1);
    b[1]=a[0];
    FOR(i,2,n+1)b[i]=a[i-1]*inv[i];
    return b;
}

inline V<mi> poly_ln(const V<mi> &a,int g){ // b=ln(a) verified with lg4725
    int n=a.size();
    assert(n),assert(a[0].val==1);
    V<mi>b=poly_conv_add(poly_diff(a),poly_inv(a,g),g);
    b.resize(n);
    return poly_intg(b);
}

inline V<mi> poly_exp(const V<mi> &a,int g){ // b=exp(a) verified with lg4726
    int n=a.size();
    assert(n);
    V<mi>b{1};
    if(a[0].val){
        mi e=0,ifac=mod-1;
        Rep(i,mod)e+=ifac,ifac*=i;
        b[0]=e^a[0].val; // check that a[0] isnt modulo
    }
    while(b.size()<a.size()){
        int m=min(b.size()<<1,a.size());
        b.resize(m);
        V<mi>c=poly_ln(b,g);
        For(i,m)c[i]=a[i]-c[i];
        ++c[0];
        b=poly_conv_add(b,c,g);
        b.resize(m);
    }
    return b;
}

inline V<mi> poly_series(const V<mi> &a,mi b0,int g){ // b[i]=¡Æ(b[j]*a[i-j]) for j>0 verified with lg4721
    assert(a.size());
    V<mi>b=a;
    b[0]=1;
    FOR(i,1,b.size())b[i]=-b[i];
    b=poly_inv(b,g);
    if(b0.val!=1)for(mi &i:b)i*=b0;
    return b;
}

inline V<mi> poly_pow(const V<mi> &_a,mi b,int g){ // c=a^(b%mod) verified with lg5245
    int n=_a.size();
    assert(n);
    V<mi>a(n);
    if(!b){
        a[0]=1;
        return a;
    }
    int i=0;
    while(i<n&&!_a[i])++i;
    if(i==n)return a;
    ll z=1ll*b.val*i;
    if(z>=n)return a;
    assert(_a[i].val==1);
    a=poly_ln(V<mi>(_a.begin()+i,_a.end()),g);
    for(mi &j:a)j*=b;
    a=poly_exp(a,g);
    V<mi>ret(z);
    ret.insert(ret.end(),a.begin(),a.begin()+n-z);
    return ret;
}

inline V<mi> poly_pow(const V<mi> &_a,ll b,int g){ // c=a^b verified with Library Checker
    int n=_a.size();
    assert(n);
    V<mi>a(n);
    if(!b){
        a[0]=1;
        return a;
    }
    int i=0;
    while(i<n&&!_a[i])++i;
    if(i==n||__int128(b)*i>=n)return a;
    a=V<mi>(_a.begin()+i,_a.end());
    mi coef=a[0],inv=1/coef;
    for(mi &j:a)j*=inv;
    a=poly_ln(a,g);
    mi _b=b%mod;
    for(mi &j:a)j*=_b;
    a=poly_exp(a,g);
    coef^=b%(mod-1);
    for(mi &j:a)j*=coef;
    ll z=b*i;
    V<mi>ret(z);
    ret.insert(ret.end(),a.begin(),a.begin()+n-z);
    return ret;
}

inline V<mi> poly_multi_pt(const V<mi> &_a,const V<mi> &b,int g){ // c[i]=a(b[i]) verified with lg5050
    assert(_a.size());
    if(b.empty())return {};
    int n=max(_a.size(),b.size());
    V<V<mi>>t(n<<2);
    auto build=[&](auto &&self,int p,int l,int r)->void{
        if(l==r){
            t[p]={1,l<b.size()?-b[r]:0};
            return;
        }
        int mid=l+r>>1;
        self(self,p<<1,l,mid);
        self(self,p<<1|1,mid+1,r);
        t[p]=poly_conv_add(t[p<<1],t[p<<1|1],g);
    };
    build(build,1,0,n-1);
    V<mi>ret(b.size());
    auto push_down=[&](auto &&self,int p,int l,int r,V<mi> c)->void{
        if(l>=b.size())return;
        if(l==r){
            ret[l]=c[0];
            return;
        }
        c.resize(r-l+1);
        int mid=l+r>>1;
        self(self,p<<1,l,mid,poly_conv_sub(c,t[p<<1|1],g));
        self(self,p<<1|1,mid+1,r,poly_conv_sub(c,t[p<<1],g));
    };
    V<mi>a=_a;
    a.resize(n+1);
    push_down(push_down,1,0,n-1,poly_conv_sub(a,poly_inv(t[1],g),g));
    return ret;
}

inline V<mi> poly_prod(const V<V<mi>> &a,int g){ // b=¡Ç(a[i])
    assert(a.size());
    auto cmp=[&](const V<mi> &x,const V<mi> &y){return x.size()>y.size();};
    priority_queue<V<mi>,V<V<mi>>,decltype(cmp)>q(cmp);
    for(const auto &i:a)q.push(i);
    while(q.size()>1){
        V<mi>x=q.top();q.pop();
        V<mi>y=q.top();q.pop();
        q.push(poly_conv_add(x,y,g));
    }
    return q.top();
}
            
inline V<mi> poly_multi_pt_sum(const V<mi> &a,int m,int g){ // b[i]=sum(a[j]^i) for i in [0,m]
    int n=a.size();
    assert(n);
    V<V<mi>>b(max(n,m));
    For(i,max(n,m))b[i]={1,-a[i]};
    V<mi>c=poly_ln(poly_prod(b,g),g);
    c.resize(m+1);
    c[0]=n;
    FOR(i,1,m+1)c[i]*=mod-i;
    return c;
}

inline pair<V<mi>,V<mi>> poly_recur(V<mi> a,V<mi> c,int g){ // build P(x)/Q(x) from sum(a^i*c^(k-i))=0 for i in [0,k] verified by lg4723
    int k=a.size();
    assert(k),assert(k+1==c.size());
    assert(c[0].val);
    a=poly_conv_add(a,c,g);
    a.resize(k);
    return {a,c};
}
inline mi poly_coef(ll m,V<mi> P,V<mi> Q,int g){ // [x^m] P(x)/Q(x) verified by abc436g
    assert(P.size()&&Q.size()),assert(Q[0].val);
    for(;m;m>>=1){
        V<mi>R=Q;
        for(int i=1;i<Q.size();i+=2)R[i]=-R[i];
        P=poly_conv_add(P,R,g),Q=poly_conv_add(Q,R,g);
        int i;
        for(i=m&1;i<P.size();i+=2)P[i>>1]=P[i];
        P.resize(i>>1);
        if(P.empty())return 0;
        for(i=0;i<Q.size();i+=2)Q[i>>1]=Q[i];
        Q.resize(i>>1);
    }
    return P[0]/Q[0];
}

inline V<mi> bernoulli(int n,int g){ // a[i]=Bernoulli[i]/i!
    assert(n>0);
    V<mi>a(n+1);
    a[1]=1;
    a=poly_exp(a,g);
    a.erase(a.begin());
    return poly_inv(a,g);
}

inline mi poly_intv_sum(V<mi> a,ll l,ll r,int g){ // ¡Æ(a[i]*x^i) for l<=x<=r
    if(l>r)return 0;
    int n=a.size();
    assert(n);
    mi fac=1;
    a.insert(a.begin(),0),++n;
    FOR(i,1,n)a[i]*=fac,fac*=i;
    a=poly_conv_sub(a,bernoulli(n,g),g);
    fac=1/fac;
    Rep(i,n)a[i]*=fac,fac*=i;
    l%=mod,r=(r+1)%mod;
    if(l<0)l+=mod;
    if(r<0)r+=mod;
    mi pl=1,pr=1,ret=0;
    For(i,n)ret+=a[i]*(pr-pl),pl*=l,pr*=r;
    return ret;
}

inline mi poly_intv_sum(V<mi> a,ll l,ll r){ // ¡Æ(a[i]*x^i) for l<=x<=r
    if(l>r)return 0;
    int n=a.size();
    assert(n);
    V<mi>S{1};
    S.reserve(n-1);
    FOR(i,1,n){
        Rep(j,i)a[j]+=(S[j]=(j?S[j-1]:0)+j*S[j])*a[i];
        S.pb(1);
    }
    l=(l-1)%mod,r%=mod;
    if(l<0)l+=mod;
    if(r<0)r+=mod;
    static V<mi>inv{0,1};
    while(inv.size()<=n)inv.pb((mod-mod/inv.size())*inv[mod%inv.size()]);
    mi cl=1,cr=1,ret=0;
    For(i,n)cl*=l+1-i,cr*=r+1-i,ret+=a[i]*(cr-cl)*inv[i+1];
    return ret;
}
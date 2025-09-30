inline V<mi> poly_conv_add(const V<mi> &_a,const V<mi> &_b,int g){ // c[k]=¡Æ(a[i]*b[j]) for i+j=k verified with lg3803
    assert(_a.size()&&_b.size());
    if(max(_a.size(),_b.size())<17){
        V<mi>c(_a.size()+_b.size()-1);
        For(i,_a.size())For(j,_b.size())c[i+j]+=_a[i]*_b[j];
        return c;
    }
    int lg=0,n=1;
    while(n<_a.size()+_b.size()-1)++lg,n<<=1;
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    static V<V<int>>btf;
    while(btf.size()<=lg){
        int n=1<<btf.size();
        btf.pb({});
        V<int>&bf=btf.back();
        bf.resize(n);
        For(i,n)bf[i]=(bf[i>>1]>>1)|((i&1)?n>>1:0);
    }
    const V<int>&bf=btf[lg];
    auto NTT=[&](V<mi> &f,mi coef){
        For(i,n)if(i<bf[i])swap(f[i],f[bf[i]]);
        for(int k=1,l=2;k<n;k<<=1,l<<=1){
            mi wn=coef^((mod-1)/l);
            for(int i=0;i<n;i+=l){
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
    For(i,n)a[i]*=b[i];
    NTT(a,mi(1)/g);
    a.resize(_a.size()+_b.size()-1);
    mi invn=mi(1)/n;
    for(mi &i:a)i*=invn;
    return a;
}

inline V<mi> poly_conv_sub(const V<mi> &_a,const V<mi> &_b,int g){ // c[k]=¡Æ(a[i]*b[j]) for i-j=k verified with gym105386H
    assert(_a.size()&&_b.size());
    V<mi>b=_b;
    reverse(ALL(b));
    b=poly_conv_add(_a,b,g);
    // (-b.size(),a.size()) -> [0,a.size())
    b.erase(b.begin(),b.begin()+_b.size()-1);
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
inline V<mi> poly_conv_mul(const V<mi> &_a,const V<mi> &_b,int g,int p,int pg=-1){ // c[k]=¡Æ(a[i]*b[j]) for i*j%p=k verified by qoj9247
    assert(_a.size()&&_b.size());
    if(!~pg)pg=find_g(p);
    assert(~pg);
    V<int>exp(p-1),lg(p);
    lg[0]=-1;
    for(int i=1,j=0;j<p-1;i=1ll*i*pg%p,++j)exp[j]=i,lg[i]=j;
    V<mi>a(p-1),b(p-1);
    FOR(i,1,_a.size())a[lg[i]]=_a[i];
    FOR(i,1,_b.size())b[lg[i]]=_b[i];
    V<mi>c=poly_conv_add(a,b,g);
    FOR(i,p-1,c.size())c[i-(p-1)]+=c[i];
    V<mi>d(p);
    d[0]=_a[0]*reduce(ALL(_b))+_b[0]*reduce(ALL(_a))-_a[0]*_b[0];
    For(i,p-1)d[exp[i]]=c[i];
    return d;
}

inline V<mi> poly_conv_div(const V<mi> &_a,const V<mi> &_b,int g,int p,int pg=-1){ // c[k]=¡Æ(a[i]*b[j]) for i/j%p=k not verified
    assert(_a.size()&&_b.size()),assert(!_b[0].val);
    V<int>inv(p);
    inv[1]=1;
    FOR(i,1,p)inv[i]=1ll*(p-p/i)*inv[p%i]%mod;
    V<mi>b(p);
    FOR(i,1,_b.size())b[inv[i]]=_b[i];
    return poly_conv_mul(_a,b,g,p,pg);
}

inline V<mi> poly_conv_and(const V<mi> &_a,const V<mi> &_b){ // c[k]=¡Æ(a[i]*b[j]) for i&j=k verified with lg4717
    assert(_a.size()&&_b.size());
    int n=1;
    while(n<max(_a.size(),_b.size()))n<<=1;
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<n;k<<=1,l<<=1)for(int i=0;i<n;i+=l)For(j,k)f[i|j]+=f[i|j|k]*coef;
    };
    FWT(a,1),FWT(b,1);
    For(i,n)a[i]*=b[i];
    FWT(a,mod-1);
    return a;
}

inline V<mi> poly_conv_or(const V<mi> &_a,const V<mi> &_b){ // c[k]=¡Æ(a[i]*b[j]) for i|j=k verified with lg4717
    assert(_a.size()&&_b.size());
    int n=1;
    while(n<max(_a.size(),_b.size()))n<<=1;
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<n;k<<=1,l<<=1)for(int i=0;i<n;i+=l)For(j,k)f[i|j|k]+=f[i|j]*coef;
    };
    FWT(a,1),FWT(b,1);
    For(i,n)a[i]*=b[i];
    FWT(a,mod-1);
    return a;
}

inline V<mi> poly_conv_xor(const V<mi> &_a,const V<mi> &_b){ // c[k]=¡Æ(a[i]*b[j]) for i^j=k verified with lg4717
    assert(_a.size()&&_b.size());
    int n=1;
    while(n<max(_a.size(),_b.size()))n<<=1;
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    auto FWT=[&](V<mi> &f,int coef){
        for(int k=1,l=2;k<n;k<<=1,l<<=1)for(int i=0;i<n;i+=l)For(j,k){
            mi x=f[i|j],y=f[i|j|k];
            f[i|j]=(x+y)*coef,f[i|j|k]=(x-y)*coef;
        }
    };
    FWT(a,1),FWT(b,1);
    For(i,n)a[i]*=b[i];
    FWT(a,mod+1>>1);
    return a;
}



inline V<mi> poly_conv_gcd(const V<mi> &_a,const V<mi> &_b){ // c[k]=¡Æ(a[i]*b[j]) for gcd(i,j)=k verified with lc418t4
    assert(_a.size()&&_b.size());
    int n=max(_a.size(),_b.size());
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    V<int>pri;
    V<bool>vis(n);
    FOR(i,2,n)if(!vis[i]){
        pri.pb(i);
        for(int k=(n-1)/i,j=k*i;k;j-=i,--k)a[k]+=a[j],b[k]+=b[j],vis[j]=true;
    }
    FOR(i,1,n)a[i]*=b[i];
    for(int i:pri)for(int j=i,k=1;j<n;j+=i,++k)a[k]-=a[j];
    a[0]=_a[0]*_b[0];
    FOR(i,1,n)a[i]+=_a[0]*_b[i]+_b[0]*_a[i];
    return a;
}

inline V<mi> poly_conv_lcm(const V<mi> &_a,const V<mi> &_b){ // c[k]=¡Æ(a[i]*b[j]) for lcm(i,j)=k not verified
    assert(_a.size()&&_b.size());
    int n=max(_a.size(),_b.size());
    V<mi>a=_a,b=_b;
    a.resize(n),b.resize(n);
    V<int>pri;
    V<bool>vis(n);
    FOR(i,2,n)if(!vis[i]){
        pri.pb(i);
        for(int j=i,k=1;j<n;j+=i,++k)a[j]+=a[k],b[j]+=b[k],vis[j]=true;
    }
    FOR(i,1,n)a[i]*=b[i];
    for(int i:pri)for(int k=(n-1)/i,j=k*i;k;j-=i,--k)a[j]-=a[k];
    a[0]=_a[0]*_b[0];
    FOR(i,1,n)a[i]+=_a[0]*_b[i]+_b[0]*_a[i];
    return a;
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
    V<mi>b(n+1),inv(n+1);
    b[1]=a[0],inv[1]=1;
    FOR(i,2,n)b[i]=a[i-1]*(inv[i]=(mod-mod/i)*inv[mod%i]);
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
    auto poly_conv_sub=[&](const V<mi> &_a,const V<mi> &_b,int g){
        assert(_b.size()),assert(_a.size()>=_b.size());
        V<mi>b=_b;
        reverse(ALL(b));
        b=poly_conv_add(_a,b,g);
        return V<mi>(b.begin()+_b.size()-1,b.end());
    };
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
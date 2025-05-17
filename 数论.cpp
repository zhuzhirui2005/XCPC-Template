// assumed that [mod<=INT_MAX] is true

template<class T>
T exgcd(const T &a,const T &b,T &x,T &y){
	if(!b){x=1,y=0;return a;}
	T g=exgcd(b,a%b,y,x);
	y-=a/b*x;
	return g;
};
template<class T>
inline T inv_exgcd(T n,T p=mod){
	// n*inv = 1 (mod p)
	// n*inv + p*k = 1
	// a*x + b*y = 1
	T inv=0,tmp=0;
	exgcd(n,p,inv,tmp);
	return inv<0?inv+p:inv;
}
template<class T>
inline ll exCRT(const V<T> &a,const V<T> &m){
    int n=a.size();
    assert(n==m.size());
    For(i,n)assert(0<=a[i]&&0<m[i]);
    function<ll(ll,ll,ll)>mul=[&](ll x,ll y,ll p=mod){
        ll z=0;
        auto add=[&](ll x,ll y){return x+y>=p?x+y-p:x+y;};
        for(x%=p;y;x=add(x,x),y>>=1)(y&1)&&(z=add(z,x));
        return z;
    };
    ll md=m[0],ret=a[0],x,y;
    FOR(i,1,n){
        ll g=exgcd(md,(ll)m[i],x,y),res=a[i]-ret%m[i];
        if(res<0)res+=m[i];
        if(res%g)return -1;
        ll mg=m[i]/g;
        if(x<0)x+=m[i];
        ret+=(x=mul(x,res/g,mg))*md;
		ret%=(md*=mg);
		if(ret<0)ret+=md;
	}
	return ret;
}

inline V<int> inverse(int n,int p=mod){
	V<int>inv(n+1);
	inv[1]=1;
	FOR(i,2,n+1)inv[i]=1ll*(p-p/i)*inv[p%i]%p;
	return inv;
}

inline V<V<int>> comb(int n,int m=-1,int p=mod){
	if(m==-1)m=n;
	if(n<m||m<0)return V<V<int>>();
    V<V<int>>C(n+1,V<int>(m+1));
    For(i,n+1){
        C[i][0]=1;
        FOR(j,1,min(i+1,m+1)){
            C[i][j]=C[i-1][j-1]+C[i-1][j];
            if(C[i][j]>=p)C[i][j]-=p;
        }
    }
    return C;
}

struct comb_table{
	int n;
	V<mi>fac,ifac;
	inline comb_table(int n_=0){n=n_,init();}
	inline void init(){
		V<mi>(n+1).swap(fac),V<mi>(n+1).swap(ifac);
		fac[0]=1;
		FOR(i,1,n+1)fac[i]=fac[i-1]*i;
		ifac[n]=1/fac[n];
		Rep(i,n)ifac[i]=ifac[i+1]*(i+1);
	}
	inline mi C(int x,int y){return x<y||y<0?0:fac[x]*ifac[y]*ifac[x-y];}
};

struct pri_table{
	int n;
	// fac[i] is the minimum prime factor of i
	V<int>fac,pri;
	inline pri_table(int n_=0){n=n_,init();}
	inline void init(){
        if(n<1)return;
		V<int>(n+1).swap(fac),V<int>().swap(pri);
		fac[1]=1;
		FOR(i,2,n+1){
			if(!fac[i])fac[i]=i,pri.pb(i);
			for(int j:pri){
				if(i*j>n)break;
				fac[i*j]=j;
				if(i%j==0)break;
			}
		}
	}
	inline bool isp(int k){return k<2?false:(fac[k]==k);}
	inline V<int> div(int k){
		assert(k<=n);
		if(k<2)return V<int>();
		V<int>ret;
		while(k>1){
			int f=fac[k];
			do k/=f;while(k%f==0);
			ret.pb(f);
		}
		return ret;
	}
};

struct mu_table{
	int n;
	V<int>mu,pri;
	V<bool>vis;
	inline mu_table(int n_=0){n=n_,init();}
	inline void init(){
        if(n<1)return;
		V<int>(n+1).swap(mu),V<int>().swap(pri),V<bool>(n+1).swap(vis);
		mu[1]=1;
		FOR(i,2,n+1){
			if(!vis[i])mu[i]=-1,pri.pb(i);
			for(int j:pri){
				if(i*j>n)break;
				vis[i*j]=true;
				if(i%j==0)break;
				mu[i*j]=-mu[i];
			}
		}
	}
	inline int get(int k){return k<1?0:mu[k];}
};

struct phi_table{
	int n;
	V<int>phi,pri;
	inline phi_table(int n_=0){n=n_,init();}
	inline void init(){
        if(n<1)return;
		V<int>(n+1).swap(phi),V<int>().swap(pri);
		phi[1]=1;
		FOR(i,2,n+1){
			if(!phi[i])phi[i]=i-1,pri.pb(i);
			for(int j:pri){
				if(i*j>n)break;
				if(i%j==0){
					phi[i*j]=phi[i]*j;
					break;
				}
				phi[i*j]=phi[i]*(j-1);
			}
		}
	}
	inline int get(int k){return k<1?0:phi[k];}
};

struct d_table{
	int n;
	V<int>cnt,d,pri;
	V<bool>vis;
	inline d_table(int n_=0){n=n_,init();}
	inline void init(){
        if(n<1)return;
		V<int>(n+1).swap(cnt),V<int>(n+1).swap(d),V<int>().swap(pri),V<bool>(n+1).swap(vis);
		cnt[1]=d[1]=1;
		FOR(i,2,n+1){
			if(!vis[i])cnt[i]=1,d[i]=2,pri.pb(i);
			for(int j:pri){
				if(i*j>n)break;
				vis[i*j]=true;
				if(i%j==0){
					int &x=cnt[i*j];
					x=cnt[i]+1;
					d[i*j]=d[i]/x*(x+1);
					break;
				}
				cnt[i*j]=1,d[i*j]=d[i]<<1;
			}
		}
	}
	inline int get(int k){return k<1?0:d[k];}
};

inline mi lagrange(int l,const V<mi> &y,int x){
    assert(y.size());
    int n=y.size();
    if(n==1)return y[0];
    if(l<=x&&x<l+n)return y[x-l];
    int r=l+n-1;
    r%=mod;if(r<0)r+=mod;
    x%=mod;if(x<0)x+=mod;
    if(r>=x&&x>r-n)return y[n-(r-x)-1];
    V<mi>ifac(n);
    ifac[0]=ifac[1]=1;
    FOR(i,2,n)ifac[i]=(mod-mod/i)*ifac[mod%i];
    FOR(i,2,n)ifac[i]*=ifac[i-1];
    V<mi>suf(n);
    suf[n-1]=1;
    REP(i,1,n)suf[i-1]=suf[i]*(x+mod-r+n-1-i);
    mi pre=1,ret=0;
    For(i,n){
        if((n-i)&1)ret+=y[i]*pre*suf[i]*ifac[i]*ifac[n-1-i];
        else ret-=y[i]*pre*suf[i]*ifac[i]*ifac[n-1-i];
        pre*=x+mod-r+n-1-i;
    }
    return ret;
}
inline mi sumexp(int n,int k){
    assert(min(n,k)>=0);
    V<int>pri;
    V<mi>pw(k+2);
    pw[0]=!k,pw[1]=1;
    FOR(i,2,k+2){
        if(!pw[i])pri.pb(i),pw[i]=mi(i)^k;
        for(int j:pri){
            if(i*j>k+1)break;
            pw[i*j]=pw[i]*pw[j];
            if(i%j==0)break;
        }
    }
    FOR(i,2-!k,k+2)pw[i]+=pw[i-1];
    return lagrange(0,pw,n);
}

/*
pre_f=sum(mu) pre_g=n pre_fg=1
pre_f=sum(phi) pre_g=n pre_fg=n*(n+1)/2
pre_f=sum(phi*id) pre_g=n*(n+1)/2 pre_fg=n*(n+1)*(2n+1)/6
*/
template<class T,class container>
T du_sieve(T n,const V<T> &pre_f,const function<T(T)> &pre_g,const function<T(T)> &pre_fg,container &h){
	if(n<pre_f.size())return pre_f[n];
	auto it=h.emplace(n,0);
	T &x=it.fi->se;
	if(it.se){
		T pre=pre_g(1);
		x=pre_fg(n);
		for(T i=2;i<=n;++i){
			T div=n/i,j=n/div,cur=pre_g(j);
			x-=(cur-pre)*du_sieve(div,pre_f,pre_g,pre_fg,h);
			i=j,pre=cur;
		}
	}
	reurn x;
}

struct vote_1{
	pii v;
	inline vote_1(){v={-1,0};}
	inline vote_1(int id,int cnt=1){v={id,cnt};}
	inline vote_1 operator+(const vote_1 &rhs){
		vote_1 ret=*this;
		if(!~ret.v.fi)ret=rhs;
		else if(~rhs.v.fi){
			if(ret.v.fi==rhs.v.fi)ret.v.se+=rhs.v.se;
			else{
				if(ret.v.se<rhs.v.se)ret={rhs.v.fi,rhs.v.se-ret.v.se};
				else ret.v.se-=rhs.v.se;
			}
		}
		return ret;
	}
};

template<int(*n)()>
struct vote{
	V<pii>v;
	inline vote(){V<pii>(n(),{-1,0}).swap(v);}
	inline vote(int id,int cnt=1){V<pii>(n(),{-1,0}).swap(v),v[0]={id,cnt};}
	inline vote operator+(const vote<n> &rhs){
		vote<n>ret=*this;
		for(pii i:rhs.v){
			if(!~i.fi)break;
			for(pii &j:ret.v)if(!~j.fi||i.fi==j.fi){
				j.fi=i.fi,j.se+=i.se;
				goto skip;
			}
			for(pii &j:ret.v)if(i.se>j.se)swap(i,j);
			for(pii &j:ret.v)j.se-=i.se;
			skip:;
		}
		return ret;
	}
};

template<int w2>
struct fp2{
    mi a,b;
    inline fp2(mi _a=0,mi _b=0):a(_a),b(_b){}
    inline fp2 operator+(mi rhs)const{return fp2(a+rhs,b);}
    inline fp2 operator-(mi rhs)const{return fp2(a-rhs,b);}
    inline fp2 operator*(mi rhs)const{return fp2(a*rhs,b*rhs);}
    inline fp2 operator/(mi rhs)const{mi inv=1/rhs;return fp2(a*inv,b*inv);}
    inline fp2 operator^(int k)const{fp2 pw=*this,ret(1);for(;k;k>>=1,pw=pw*pw)if(k&1)ret=ret*pw;return ret;}
    inline fp2& operator+=(mi rhs){a+=rhs;return *this;}
    inline fp2& operator-=(mi rhs){a-=rhs;return *this;}
    inline fp2& operator*=(mi rhs){a*=rhs,b*=rhs;return *this;}
    inline fp2& operator/=(mi rhs){mi inv=1/rhs;a*=inv,b*=inv;return *this;}
    inline fp2& operator^=(int k){fp2 tmp(1),base=*this;for(;k;k>>=1,base*=base)if(k&1)tmp*=base;return *this=tmp;}
    inline fp2 operator+(const fp2&rhs)const{return fp2(a+rhs.a,b+rhs.b);}
    inline fp2 operator-(const fp2&rhs)const{return fp2(a-rhs.a,b-rhs.b);}
    inline fp2 operator*(const fp2&rhs)const{return fp2(a*rhs.a+b*rhs.b*w2,a*rhs.b+rhs.a*b);}
    inline fp2 operator/(const fp2&rhs)const{assert(rhs.a.val||rhs.b.val);mi inv=1/(rhs.a*rhs.a-rhs.b*rhs.b*w2);return fp2((a*rhs.a-b*rhs.b*w2)*inv,(rhs.a*b-a*rhs.b)*inv);}
    inline fp2& operator+=(const fp2&rhs){a+=rhs.a,b+=rhs.b;return *this;}
    inline fp2& operator-=(const fp2&rhs){a-=rhs.a,b-=rhs.b;return *this;}
    inline fp2& operator*=(const fp2&rhs){mi x=a*rhs.a+b*rhs.b*w2,y=a*rhs.b+rhs.a*b;a=x,b=y;return *this;}
    inline fp2& operator/=(const fp2&rhs){assert(rhs.a.val||rhs.b.val);mi inv=1/(rhs.a*rhs.a-rhs.b*rhs.b*w2);mi x=(a*rhs.a-b*rhs.b*w2)*inv,y=(rhs.a*b-a*rhs.b)*inv;a=x,b=y;return *this;}
    inline fp2 operator-()const{return fp2(-a,-b);}
    friend fp2 operator+(mi lhs,const fp2&rhs){return fp2(lhs+rhs.a,rhs.b);}
    friend fp2 operator-(mi lhs,const fp2&rhs){return fp2(lhs-rhs.a,-rhs.b);}
    friend fp2 operator*(mi lhs,const fp2&rhs){return fp2(lhs*rhs.a,lhs*rhs.b);}
    friend fp2 operator/(mi lhs,const fp2&rhs){assert(rhs.a.val||rhs.b.val);mi inv=1/(rhs.a*rhs.a-rhs.b*rhs.b*w2);return fp2(lhs*rhs.a*inv,-lhs*rhs.b*inv);}
};
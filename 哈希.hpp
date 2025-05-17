template<int base=2333>
struct mhsh{
	// 0-indexed
	V<ull>bs,h;
	inline mhsh(){}
	inline mhsh(const string &s){
		bs.reserve(s.size()),h.reserve(s.size());
		bs.pb(1),h.pb(s[0]);
		FOR(i,1,s.size())bs.pb(bs.back()*base),h.pb(h.back()*base+s[i]);
	}
	inline mhsh(const V<int> &v){
        bs.reserve(v.size()),h.reserve(v.size());
        bs.pb(1),h.pb(v[0]);
        FOR(i,1,v.size())bs.pb(bs.back()*base),h.pb(h.back()*base+v[i]);
    }
	inline ull get(int l,int r){
		assert(0<=l),assert(l<=r),assert(r<h.size());
		return h[r]-(l?h[l-1]*bs[r-l+1]:0);
	}
	inline int lcp(int x,int y){
		assert(0<=min(x,y)),assert(max(x,y)<h.size());
		int l=1,r=h.size()-max(x,y),ret=0;
		while(l<=r){
			int mid=l+r>>1;
			if(get(x,x+mid-1)==get(y,y+mid-1))l=mid+1,ret=mid;
			else r=mid-1;
		}
		return ret;
	}
};

template<int base=2337,int mod=998244853>
struct modhsh{
    // 0-indexed
    V<ull>bs,h;
    inline modhsh(){}
    inline modhsh(const string &s){
        bs.reserve(s.size()),h.reserve(s.size());
        bs.pb(1),h.pb(s[0]);
        FOR(i,1,s.size())bs.pb(bs.back()*base%mod),h.pb((h.back()*base+s[i])%mod);
    }
    inline modhsh(const V<int> &v){
        bs.reserve(v.size()),h.reserve(v.size());
        bs.pb(1),h.pb(v[0]);
        FOR(i,1,v.size())bs.pb(bs.back()*base%mod),h.pb((h.back()*base+v[i])%mod);
    }
    inline ull get(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<h.size());
        ull ret=h[r]+mod-(l?h[l-1]*bs[r-l+1]%mod:0);
        return ret>=mod?ret-mod:ret;
    }
};
template<int base1=2337,int mod1=998244853,int base2=2333,int mod2=1'000'000'009>
struct dmhsh{
    // 0-indexed
    modhsh<base1,mod1>hsh1;
    modhsh<base2,mod2>hsh2;
    inline dmhsh(const string &s){
        hsh1=modhsh<base1,mod1>(s),hsh2=modhsh<base2,mod2>(s);
    }
    inline dmhsh(const V<int> &v){
        hsh1=modhsh<base1,mod1>(v),hsh2=modhsh<base2,mod2>(v);
    }
    inline pair<ull,ull> get(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<hsh1.h.size());
        return {hsh1.get(l,r),hsh2.get(l,r)};
    }
    inline int lcp(int x,int y){
		assert(0<=min(x,y)),assert(max(x,y)<hsh2.h.size());
		int l=1,r=hsh2.h.size()-max(x,y),ret=0;
		while(l<=r){
			int mid=l+r>>1;
			if(get(x,x+mid-1)==get(y,y+mid-1))l=mid+1,ret=mid;
			else r=mid-1;
		}
		return ret;
	}
};

mt19937 rnd(time(0));
inline int genPri(int l,int r){
    auto isp=[&](int k){
        if(k<2)return false;
        for(int i=2;i*i<=k;++i)if(k%i==0)return false;
        return true;
    };
    int p=uniform_int_distribution<int>(l,r)(rnd);
    while(!isp(p))++p;
    return p;
};
struct rndhsh{
    // 0-indexed
    int base,mod;
    V<ull>bs,h;
    inline rndhsh(){base=genPri(2,1e5),mod=genPri(2,1e9);}
    inline rndhsh(const string &s){
        bs.reserve(s.size()),h.reserve(s.size());
        bs.pb(1),h.pb(s[0]);
        FOR(i,1,s.size())bs.pb(bs.back()*base%mod),h.pb((h.back()*base+s[i])%mod);
    }
    inline rndhsh(const V<int> &v){
        bs.reserve(v.size()),h.reserve(v.size());
        bs.pb(1),h.pb(v[0]);
        FOR(i,1,v.size())bs.pb(bs.back()*base%mod),h.pb((h.back()*base+v[i])%mod);
    }
    inline ull get(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<h.size());
        ull ret=h[r]+mod-(l?h[l-1]*bs[r-l+1]%mod:0);
        return ret>=mod?ret-mod:ret;
    }
};
struct drhsh{
    // 0-indexed
    rndhsh hsh1;
    rndhsh hsh2;
    inline drhsh(const string &s){
        hsh1=rndhsh(s),hsh2=rndhsh(s);
    }
    inline drhsh(const V<int> &v){
        hsh1=rndhsh(v),hsh2=rndhsh(v);
    }
    inline pair<ull,ull> get(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<hsh1.h.size());
        return {hsh1.get(l,r),hsh2.get(l,r)};
    }
    inline int lcp(int x,int y){
		assert(0<=min(x,y)),assert(max(x,y)<hsh2.h.size());
		int l=1,r=hsh2.h.size()-max(x,y),ret=0;
		while(l<=r){
			int mid=l+r>>1;
			if(get(x,x+mid-1)==get(y,y+mid-1))l=mid+1,ret=mid;
			else r=mid-1;
		}
		return ret;
	}
};
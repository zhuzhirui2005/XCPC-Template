inline int lcs(const string &a,const string &b){
	if(a.empty()||b.empty())return 0;
	int n=a.size(),m=b.size(),k=(n+62)/63;
	V<ull>f(k);
	char mn=*min_element(ALL(a)),mx=*max_element(ALL(a));
	V<V<ull>>g(mx-mn+1,V<ull>(k));
	For(i,n)g[a[i]-mn][i/63]|=1ull<<i%63;
	for(char i:b){
		if(i<mn||i>mx)continue;
        i-=mn;
		ull z=1;
		For(j,k){
			ull x=f[j],y=f[j]|g[i][j];
			((x<<=1)|=z)+=(~y)&((1ull<<63)-1);
			f[j]=x&y,z=x>>63;
		}
	}
	return accumulate(ALL(f),0,[&](int x,ull y){return x+__builtin_popcountll(y);});
}
template<class T>
inline int lcs(const V<T> &a,const V<T> &b){
	if(a.empty()||b.empty())return 0;
	int n=a.size(),m=b.size(),k=(n+62)/63;
	disc<T>d(a);
	V<ull>f(k);
	V<V<ull>>g(d.size(),V<ull>(k));
	For(i,n)g[d.query(a[i])][i/63]|=1ull<<i%63;
	for(const T &i:b){
		auto it=lower_bound(ALL(d.d),i);
		if(it==d.d.end()||*it!=i)continue;
		i=it-d.d.begin();
		ull z=1;
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
	inline subseq_table(const string &v){
		int n=v.size();
		V<V<int>>(128).swap(nxt);
		For(i,n){
			assert(v[i]>=0&&v[i]<128);
			nxt[v[i]].pb(i);
		}
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

template<class T>
struct subseq_Table{
	genUID<T>g;
	V<V<int>>nxt;
	inline subseq_Table(const V<T> &v){
		int n=v.size();
		For(i,n){
			int k=g.get_id(v[i]);
			if(k>=nxt.size())nxt.pb(V<int>());
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

struct manacher{
    int n;
    V<int>p;
    inline manacher(const string &s){
        n=s.size();
        p.assign(n<<1|1,1);
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
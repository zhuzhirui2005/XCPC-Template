template<class T>
inline V<pii> cart_seq(const V<T> &v,function<bool(T,T)>cmp=[](T x,T y){return x>y;}){
	int n=v.size();
	V<pii>ret(n,pii(-1,n));
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
	int n,rt;
	V<V<int>>to;
	inline void resize(int n_){V<V<int>>(n=n_).swap(to);}
	inline lca_table(int n_=0){resize(n_);}
	inline void add_edge(int x,int y){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
		to[x].pb(y),to[y].pb(x);
	}
	inline lca_table(const V<V<int>>&to_){n=(to=to_).size();init();}
	V<int>dep,fa,siz,son,top;
	inline void init(int _rt=0){
		rt=_rt;
		V<int>(n).swap(dep),V<int>(n).swap(fa),V<int>(n).swap(siz),V<int>(n,-1).swap(son);
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
		V<int>(n,-1).swap(top);
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

struct tree_chain{
	int n,rt;
	V<V<int>>to;
	inline void resize(int n_){V<V<int>>(n=n_).swap(to);}
	inline tree_chain(int n_=0){resize(n_);}
	inline void add_edge(int x,int y){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(x!=y);
		to[x].pb(y),to[y].pb(x);
	}
	inline tree_chain(const V<V<int>>&to_){n=(to=to_).size();init();}
	V<int>dep,fa,rev,seg,siz,son,top;
	inline void init(int _rt=0){
		rt=_rt;
		V<int>(n).swap(dep),V<int>(n).swap(fa),V<int>(n).swap(siz),V<int>(n,-1).swap(son);
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
		int cnt=0;
		V<int>(n).swap(rev),V<int>(n).swap(seg),V<int>(n,-1).swap(top);
		dfs1(rt,-1);
		function<void(int,int)>dfs2=[&](int p,int k){
			seg[p]=cnt,rev[cnt++]=p,top[p]=k;
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
        return rev[seg[p]-k];
    }
	inline V<pii> path(int x,int y,bool dir=0){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n);
		V<pii>ret,ter;
		bool rv=0;
		while(top[x]!=top[y]){
			if(dep[top[x]]<dep[top[y]])rv^=1,swap(x,y);
			if(dir){
				if(rv)ter.eb(seg[top[x]],seg[x]);
				else ret.eb(seg[x],seg[top[x]]);
			}
			else (rv?ter:ret).eb(seg[top[x]],seg[x]);
			x=fa[top[x]];
		}
		if(dep[x]>dep[y])rv^=1,swap(x,y);
		if(dir){
			if(rv)ter.eb(seg[y],seg[x]);
			else ret.eb(seg[x],seg[y]);
		}
		else (rv?ret:ter).eb(seg[x],seg[y]);
		reverse(ALL(ter));
		ret.insert(ret.end(),ALL(ter));
		return ret;
	}
};
inline void virt_tree(V<int> &p,const tree_chain &tc,V<V<int>> &to){
    sort(ALL(p),[&](int x,int y){return tc.seg[x]<tc.seg[y];});
    p.erase(unique(ALL(p)),p.end());
    auto add_edge=[&](int x,int y){to[x].pb(y),to[y].pb(x);};
    V<int>st;
    for(int i:p){
        if(st.size()){
            int anc=tc.lca(i,st.back());
            if(anc!=st.back()){
                while(st.size()>1&&tc.seg[anc]<tc.seg[st[st.size()-2]])add_edge(st[st.size()-2],st.back()),st.qb();
                if(st.size()==1||tc.seg[anc]>tc.seg[st[st.size()-2]])V<int>().swap(to[anc]),add_edge(anc,st.back()),st.back()=anc;
                else add_edge(anc,st.back()),st.qb();
            }
        }
        V<int>().swap(to[i]),st.pb(i);
    }
    while(st.size()>1)add_edge(st[st.size()-2],st.back()),st.qb();
}

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
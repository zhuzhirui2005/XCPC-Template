inline V<ll> bfs01(int n,int s,const V<V<pii>> &to){
	assert(0<=n),assert(0<=s),assert(s<n),assert(to.size()<=n);
	for(const V<pii> &i:to)
        for(const pii &j:i)
            assert(0<=min(j.fi,j.se)),assert(j.fi<n);
    V<ll>dis(n,infl);
    dis[s]=0;
    deque<int>q;
    q.pb(s);
    V<bool>vis(n); // added vis to prevent an obvious error
    while(q.size()){
        int p=q.front();q.qf();
        if(vis[p])continue;
        vis[p]=true;
        for(const pii &i:to[p])if(ckmin(dis[i.fi],dis[p]+i.se))i.se?q.pb(i.fi):q.pf(i.fi);
    }
    for(ll &i:dis)if(i==infl)i=-1;
    return dis;
}

template<class T>
inline V<ll> dijkstra(int n,int s,const V<V<pair<int,T>>> &to,ll null=-1){
    V<ll>dis(n,infl);
    dis[s]=0;
    typedef pair<int,ll> pil;
    auto cmp=[&](const pil &x,const pil &y){return x.se>y.se;};
    priority_queue<pil,V<pil>,decltype(cmp)>q(cmp);
    q.emplace(s,0);
    V<bool>vis(n);
    while(q.size()){
        int p=q.top().fi;q.pop();
        if(vis[p])continue;
        vis[p]=true;
        for(const auto &[i,j]:to[p])if(ckmin(dis[i],dis[p]+j)&&!vis[i])q.emplace(i,dis[i]);
    }
    for(ll &i:dis)if(i==infl)i=null;
    return dis;
}

inline V<array<int,3>> kruskal(int n,const V<V<pii>> &to,function<bool(const array<int,3> &,const array<int,3> &)>cmp=[](const array<int,3> x,const array<int,3> &y){return x[2]<y[2];}){
	assert(0<=n),assert(to.size()<=n);
	for(const V<pii> &i:to)for(const pii &j:i)assert(j.fi<n);
	V<array<int,3>>e;
	For(i,to.size())for(const pii &j:to[i])assert(0<=j.fi),assert(j.fi<n),e.pb({i,j.fi,j.se});
	sort(ALL(e),cmp);
	dsu d(n);
	V<array<int,3>>ret;
	for(auto &i:e)if(d.merge(i[0],i[1]))ret.pb(i);
	return ret;
}

struct ring{
    int clr;
    V<int>id;
    V<V<int>>scc,to;
    inline void init(const V<V<int>>&to){
        int cnt=clr=0,n=to.size();
        V<bool>cur(n);
        V<int>dfn(n),low(n);
        V<int>(n,-1).swap(id),V<V<int>>().swap(scc);
		stack<int>st;
        function<void(int)>tarjan=[&](int p){
            cur[p]=true;
            dfn[p]=low[p]=++cnt;
            st.push(p);
            for(int i:to[p]){
                assert(0<=i&&i<n);
                if(!dfn[i])tarjan(i),ckmin(low[p],low[i]);
                else if(cur[i])ckmin(low[p],dfn[i]);
            }
            if(dfn[p]==low[p]){
                scc.pb(V<int>());
                int k;
                do{
                    k=st.top();st.pop();
                    cur[k]=false,id[k]=clr,scc[clr].pb(k);
                }while(k!=p);
                ++clr;
            }
        };
        For(i,n)if(!dfn[i])tarjan(i);
        V<int>lst(clr,-1);
	    V<V<int>>(clr).swap(this->to);
		For(i,clr){
			lst[i]=i;
			for(int j:scc[i])for(int k:to[j])if(lst[id[k]]!=i)lst[id[k]]=i,this->to[i].pb(id[k]);
		}
    }
    inline ring(const V<V<int>>&to){init(to);}
    inline ring(){}
};

struct vDCC{
    int clr;
    V<bool>cut;
    V<V<int>>dcc,to;
    inline void init(const V<V<int>>&to){
        int cnt=0,n=clr=to.size();
        V<int>dfn(n),low(n);
        V<bool>(n).swap(cut),V<V<int>>().swap(dcc);
        V<V<int>>(n).swap(this->to);
        For(i,n)
			if(!dfn[i]){
				stack<int>st;
		        function<void(int,int)>tarjan=[&](int p,int fa){
		            dfn[p]=low[p]=++cnt;
		            int flag_son=0;
		            st.push(p);
		            for(int i:to[p]){
		                assert(0<=i&&i<n);
		                if(!dfn[i]){
							tarjan(i,p),ckmin(low[p],low[i]);
							if(low[i]>=dfn[p]){
								if(fa!=-1||flag_son++)cut[p]=true;
                                this->dcc.pb(V<int>()),this->to.pb(V<int>());
                                int k;
                                do{
                                    k=st.top();st.pop();
                                    this->dcc.back().pb(k);
                                    this->to[k].pb(clr),this->to[clr].pb(k);
                                }while(k!=i);
                                this->dcc.back().pb(p);
                                this->to[p].pb(clr),this->to[clr++].pb(p);
							}
						}
		                else ckmin(low[p],dfn[i]);
		            }
                    if(!~fa&&!flag_son)this->dcc.pb({p});
		        };
				tarjan(i,-1);
			}
    }
    inline vDCC(const V<V<int>>&to){init(to);}
    inline vDCC(){}
};

struct eDCC{
    int clr;
    V<V<int>>dcc,to;
    V<int>id;
    inline void init(const V<V<int>>&to){
        int cnt=clr=0,n=to.size();
        V<int>dfn(n),low(n);
        V<V<int>>().swap(dcc),V<int>(n,-1).swap(id);
        stack<int>st;
        function<void(int,int)>tarjan=[&](int p,int fa){
            dfn[p]=low[p]=++cnt;
            bool flag=false;
            st.push(p);
            for(int i:to[p]){
            	if(i!=fa){
            		if(!dfn[i])tarjan(i,p),ckmin(low[p],low[i]);
            		else ckmin(low[p],dfn[i]);
				}
				if(i==fa){
					if(flag)ckmin(low[p],dfn[i]);
					else flag=true;
				}
			}
			if(dfn[p]<=low[p]){
				dcc.pb(V<int>());
				int k;
				do{
					k=st.top();st.pop();
					id[k]=clr,dcc[clr].pb(k);
				}while(k!=p);
				++clr;
			}
        };
        For(i,n)if(!dfn[i])tarjan(i,-1);
        V<int>lst(clr,-1);
	    V<V<int>>(clr).swap(this->to);
		For(i,clr){
			lst[i]=i;
			for(int j:dcc[i])for(int k:to[j])if(lst[id[k]]!=i)lst[id[k]]=i,this->to[i].pb(id[k]);
		}
    }
    inline eDCC(const V<V<int>>&to){init(to);}
    inline eDCC(){}
};

struct range_2sat{
	int n;
	V<V<int>>to;
	inline int idx(int l,int r){return (l+r|l!=r)>>1;}
	#define p idx(l,r)
	inline void resize(int n_){
		n=n_;
		V<V<int>>((n<<1)+(n-1<<2)).swap(to);
		function<int(int,int,int)>build_dw=[&](int l,int r,int k){
			if(l==r)return (k&1)*n+l;
	        int mid=l+r>>1;
	        to[(n<<1)+k*(n-1)+p].pb(build_dw(l,mid,k));
			to[(n<<1)+k*(n-1)+p].pb(build_dw(mid+1,r,k));
	        return (n<<1)+k*(n-1)+p;
		};
		build_dw(0,n-1,0),build_dw(0,n-1,1);
		function<int(int,int,int)>build_up=[&](int l,int r,int k){
			if(l==r)return (k&1)*n+r;
	        int mid=l+r>>1;
	        to[build_up(l,mid,k)].pb((n<<1)+k*(n-1)+p);
			to[build_up(mid+1,r,k)].pb((n<<1)+k*(n-1)+p);
			return (n<<1)+k*(n-1)+p;
		};
		build_up(0,n-1,2),build_up(0,n-1,3);
	}
	inline range_2sat(){}
	inline range_2sat(int n_){resize(n_);}
	inline V<int> range_dw(int ql,int qr,int k){
		V<int>ret;
		function<void(int,int)>dfs=[&](int l,int r){
			if(ql<=l&&r<=qr){
				if(l==r)ret.pb(k*n+l);
				else ret.pb((n<<1)+k*(n-1)+p);
				return;
			}
			int mid=l+r>>1;
			if(ql<=mid)dfs(l,mid);
			if(qr>mid)dfs(mid+1,r);
		};
		dfs(0,n-1);
		return ret;
	}
	inline V<int> range_up(int ql,int qr,int k){
		V<int>ret;
		function<void(int,int)>dfs=[&](int l,int r){
			if(ql<=l&&r<=qr){
				if(l==r)ret.pb(k*n+r);
				else ret.pb((n<<1)+(k+2)*(n-1)+p);
				return;
			}
			int mid=l+r>>1;
			if(ql<=mid)dfs(l,mid);
			if(qr>mid)dfs(mid+1,r);
		};
		dfs(0,n-1);
		return ret;
	}
	#undef p
    inline void link_pp(int x,int y,bool op_x,bool op_y,bool rev=true){
		to[op_x*n+x].pb(op_y*n+y);
		if(rev)to[(op_y^1)*n+y].pb((op_x^1)*n+x);
	}
	inline void link_pr(int x,int yl,int yr,bool op_x,bool op_y,bool rev=true){
		for(int y:range_dw(yl,yr,op_y))to[op_x*n+x].pb(y);
		if(rev)for(int y:range_up(yl,yr,op_y^1))to[y].pb((op_x^1)*n+x);
	}
	inline void link_rp(int xl,int xr,int y,bool op_x,bool op_y,bool rev=true){
		for(int x:range_up(xl,xr,op_x))to[x].pb(op_y*n+y);
		if(rev)for(int x:range_dw(xl,xr,op_x^1))to[(op_y^1)*n+y].pb(x);
	}
	inline void link_rr(int xl,int xr,int yl,int yr,bool op_x,bool op_y,bool rev=true){
		V<int>X=range_up(xl,xr,op_x);
		for(int y:range_dw(yl,yr,op_y))for(int x:X)to[x].pb(y);
		if(rev){
			V<int>Y=range_up(yl,yr,op_y^1);
			for(int x:range_dw(xl,xr,op_x^1))for(int y:Y)to[y].pb(x);
		}
	}
};
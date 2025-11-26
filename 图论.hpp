inline V<V<int>> graph(int n,const V<V<int>> &e,int dir=1,int d=0){
    assert(0<=n);
	V<V<int>>to(n);
    for(const V<int> &i:e){
        assert(i.size()==2),assert(0<=min(i[0]-d,i[1]-d)),assert(max(i[0]-d,i[1]-d)<n);
        if(dir&1)to[i[0]-d].push_back(i[1]-d);
        if(dir&2)to[i[1]-d].push_back(i[0]-d);
    }
    return to;
}

inline V<V<pii>> graphw(int n,const V<V<int>> &e,int dir=1,int d=0){
    assert(0<=n);
	V<V<pii>>to(n);
    for(const V<int> &i:e){
        assert(i.size()==3),assert(0<=min(i[0]-d,i[1]-d)),assert(max(i[0]-d,i[1]-d)<n);
        if(dir&1)to[i[0]-d].emplace_back(i[1]-d,i[2]);
        if(dir&2)to[i[1]-d].emplace_back(i[0]-d,i[2]);
    }
    return to;
}

template<class T>
inline V<V<int>> delw(const V<V<pair<int,T>>> &_to){
	V<V<int>>to(_to.size());
	For(i,_to.size()){
		to[i].resize(_to[i].size());
		For(j,_to[i].size())to[i][j]=_to[i][j].fi;
	}
    return to;
}

inline V<int> bfs(int n,int s,const V<V<int>> &to){
	assert(0<=n),assert(0<=s),assert(s<n),assert(to.size()<=n);
    V<int>dis(n,-1);
    dis[s]=0;
    queue<int>q;
    q.push(s);
    while(q.size()){
        int p=q.front();q.pop();
        for(int i:to[p])if(!~dis[i])dis[i]=dis[p]+1,q.push(p);
    }
    return dis;
}

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
inline V<V<int>> bfs2d(int n,int m,int sx,int sy,const V<V<T>> &grid,const T &can){
    V<V<int>>dis(n,V<int>(m,-1));
    dis[sx][sy]=0;
    queue<pii>q;
    q.emplace(sx,sy);
    while(q.size()){
        auto [x,y]=q.front();q.pop();
        For(i,4){
            int mx=x+dir[i][0],my=y+dir[i][1];
            if(mx>=0&&mx<n&&my>=0&&my<m&&grid[mx][my]==can&&dis[mx][my]==-1)dis[mx][my]=dis[x][y]+1,q.emplace(mx,my);
        }
    }
    return dis;
}
template<class T>
inline V<V<int>> bfs2d_multi(int n,int m,V<pii>s,const V<V<T>> &grid,const T &can){
    V<V<int>>dis(n,V<int>(m,-1));
    queue<pii>q;
    for(const auto &i:s)dis[i.fi][i.se]=0,q.push(i);
    while(q.size()){
        auto [x,y]=q.front();q.pop();
        For(i,4){
            int mx=x+dir[i][0],my=y+dir[i][1];
            if(mx>=0&&mx<n&&my>=0&&my<m&&grid[mx][my]==can&&dis[mx][my]==-1)dis[mx][my]=dis[x][y]+1,q.emplace(mx,my);
        }
    }
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

inline V<array<int,3>> kruskal(int n,V<array<int,3>> &e,function<bool(const array<int,3> &,const array<int,3> &)>cmp=[](const array<int,3> x,const array<int,3> &y){return x[2]<y[2];}){
	assert(n>=0);
	for(const auto &[u,v,w]:e)assert(0<=u),assert(u<n),assert(0<=v),assert(v<n);
	sort(ALL(e),cmp);
	dsu d(n);
	V<array<int,3>>ret;
	for(auto &i:e)if(d.merge(i[0],i[1]))ret.pb(i);
	return ret;
}

inline pair<V<V<int>>,V<int>> kruskal_tree(int n,V<array<int,3>> &e,function<bool(const array<int,3> &,const array<int,3> &)>cmp=[](const array<int,3> x,const array<int,3> &y){return x[2]<y[2];}){
    int cnt=n;
    dsu d(n+n-1);
    V<V<int>>to(n+n-1);
    V<int>val(n+n-1);
    sort(ALL(e),cmp);
    for(const auto &i:e){
        int fx=d.find(i[0]),fy=d.find(i[1]);
        if(fx!=fy){
            d.fa[fx]=d.fa[fy]=cnt;
            to[cnt].pb(fx),to[cnt].pb(fy);
            val[cnt++]=i[2];
        }
    }
    assert(cnt==n+n-1);
    return {to,val};
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

struct ring_id{
    int clr;
    V<int>id;
    inline void init(const V<V<int>>&to){
        int cnt=clr=0,n=to.size();
        V<bool>cur(n);
        V<int>dfn(n),low(n);
        V<int>(n,-1).swap(id);
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
                int k;
                do{
                    k=st.top();st.pop();
                    cur[k]=false,id[k]=clr;
                }while(k!=p);
                ++clr;
            }
        };
        For(i,n)if(!dfn[i])tarjan(i);
    }
    inline ring_id(const V<V<int>>&to){init(to);}
    inline ring_id(){}
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

inline mi matrix_tree_prod(int n,const V<array<int,3>> &to,int dir=3,int rt=0){
    // dir = 1 为外向生成树，2 为内向生成树，3 为生成树
    assert(n>0);
    for(const auto &[u,v,w]:to)assert(0<=u),assert(u<n),assert(0<=v),assert(v<n);
    if(n==1)return 1;
    matrix<mi>kh(n-1,n-1);
    for(const auto &[u,v,w]:to)if(u!=v){
        int u_=u,v_=v;
        if(u>rt)--u_;
        if(v>rt)--v_;
        if(dir&1){
            if(u!=rt&&v!=rt)kh[u_][v_]-=w;
            if(v!=rt)kh[v_][v_]+=w;
        }
        if(dir&2){
            if(v!=rt&&u!=rt)kh[v_][u_]-=w;
            if(u!=rt)kh[u_][u_]+=w;
        }
    }
    return kh.det();
}

inline mi matrix_tree_sum(int n,const V<array<int,3>> &to,int dir=3,int rt=0){
    // dir = 1 为外向生成树，2 为内向生成树，3 为生成树
    assert(n>0);
    for(const auto &[u,v,w]:to)assert(0<=u),assert(u<n),assert(0<=v),assert(v<n);
    if(n==1)return 0;
    // 对 (1 + w_i x) 在 mod x^2 意义下跑矩阵树，一次项系数就是 sum(w_i)
    struct ans_t:pair<mi,mi>{
        using pair<mi,mi>::pair;
        inline ans_t operator+(const ans_t &rhs){return {fi+rhs.fi,se+rhs.se};}
        inline ans_t operator-(const ans_t &rhs){return {fi-rhs.fi,se-rhs.se};}
        inline ans_t operator*(const ans_t &rhs){return {fi*rhs.se+rhs.fi*se,se*rhs.se};}
        inline ans_t operator/(const ans_t &rhs){mi inv=1/rhs.se;return {(fi*rhs.se-rhs.fi*se)*inv*inv,se*inv};}
    };
    V kh(n-1,V<ans_t>(n-1));
    for(const auto &[u,v,w]:to)if(u!=v){
        int u_=u,v_=v;
        if(u>rt)--u_;
        if(v>rt)--v_;
        if(dir&1){
            if(u!=rt&&v!=rt)kh[u_][v_]=kh[u_][v_]-ans_t(w,1);
            if(v!=rt)kh[v_][v_]=kh[v_][v_]+ans_t(w,1);
        }
        if(dir&2){
            if(v!=rt&&u!=rt)kh[v_][u_]=kh[v_][u_]-ans_t(w,1);
            if(u!=rt)kh[u_][u_]=kh[u_][u_]+ans_t(w,1);
        }
    }
    ans_t ret={0,1};
    For(i,n-1){
        if(!kh[i][i].se)FOR(j,i+1,n-1)if(kh[j][i].se!=0){
            ret.fi=-ret.fi,ret.se=-ret.se;
            swap(kh[i],kh[j]);
            break;
        }
        if(!kh[i][i].se)return 0;
        ans_t inv=ans_t(0,1)/kh[i][i];
        ret=ret*kh[i][i];
        FOR(j,i+1,n-1){
            ans_t coef=kh[j][i]*inv;
            FOR(k,i,n-1)kh[j][k]=kh[j][k]-coef*kh[i][k];
        }
    }
    return ret.fi;
}

inline mi matrix_tree_sum(int n,int k,const V<array<int,3>> &to,int dir=3,int rt=0){
    // dir = 1 为外向生成树，2 为内向生成树，3 为生成树
    assert(n>0),assert(k>=0);
    for(const auto &[u,v,w]:to)assert(0<=u),assert(u<n),assert(0<=v),assert(v<n);
    if(n==1)return !k;
    comb_table ct(k);
    V kh(n-1,V(n-1,V<mi>(k+1)));
    for(const auto &[u,v,w]:to)if(u!=v){
        int u_=u,v_=v;
        if(u>rt)--u_;
        if(v>rt)--v_;
        V<mi>pw(k+1);
        pw[0]=1;
        For(i,k)pw[i+1]=pw[i]*w;
        if(dir&1){
            if(u!=rt&&v!=rt)For(i,k+1)kh[u_][v_][i]-=pw[i];
            if(v!=rt)For(i,k+1)kh[v_][v_][i]+=pw[i];
        }
        if(dir&2){
            if(v!=rt&&u!=rt)For(i,k+1)kh[v_][u_][i]-=pw[i];
            if(u!=rt)For(i,k+1)kh[u_][u_][i]+=pw[i];
        }
    }
    V<mi>ret(k+1);
    ret[0]=1;
    For(i,n-1){
        if(!kh[i][i][0])FOR(j,i+1,n-1)if(kh[j][i][0]!=0){
            for(mi &l:ret)l=-l;
            swap(kh[i],kh[j]);
            break;
        }
        if(!kh[i][i][0])return 0;
        V<mi>inv(k+1);
        inv[0]=1/kh[i][i][0];
        FOR(j,1,k+1){
            FOR(l,1,j+1)inv[j]+=ct.C(j,l)*kh[i][i][l]*inv[j-l];
            inv[j]*=-inv[0];
        }
        Rep(j,k+1){
            FOR(l,1,k+1-j)ret[j+l]+=ct.C(j+l,l)*ret[j]*kh[i][i][l];
            ret[j]*=kh[i][i][0];
        }
        FOR(j,i+1,n-1){
            V<mi>coef(k+1);
            For(l,k+1)For(m,k+1-l)coef[l+m]+=ct.C(l+m,m)*kh[j][i][l]*inv[m];
            FOR(l,i,n-1)For(m,k+1)For(o,k+1-m)kh[j][l][m+o]-=ct.C(m+o,m)*coef[m]*kh[i][l][o];
        }
    }
    return ret[k];
}
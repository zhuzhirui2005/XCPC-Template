struct maxflow{
    V<int>dep;
	V<pil>e;
	V<V<int>>hd;
	int n,S,T;
	inline void add_edge(int x,int y,ll z){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(z>=0);
		hd[x].pb(e.size()),e.eb(y,z),hd[y].pb(e.size()),e.eb(x,0);
	}
	inline maxflow(int _n=0,int _S=-1,int _T=-1){
		V<V<int>>(n=_n).swap(hd);
		S=_S,T=_T;
	}
    template<class U>
	inline maxflow(const V<V<pair<int,U>>> &to,int _S=-1,int _T=-1){
        V<V<int>>(n=to.size()).swap(hd);
		For(i,n)for(const auto &[j,k]:to[i])add_edge(i,j,k);
		S=_S,T=_T;
	}
	inline ll dinic(int s=-1,int t=-1){
        if(!~s)s=S;
        if(!~t)t=T;
		assert(~s),assert(~t);
		auto bfs=[&](){
			dep.assign(n,0);
			dep[s]=1;
			queue<int>q;
			q.push(s);
			while(q.size()){
				int p=q.front();q.pop();
				for(int i:hd[p])if(e[i].se&&!dep[e[i].fi])dep[e[i].fi]=dep[p]+1,q.push(e[i].fi);
			}
			return dep[t];
		};
		function<ll(int,ll)>dfs=[&](int p,ll lim){
			if(p==t)return lim;
			ll sum=0;
			for(int i:hd[p])if(e[i].se&&dep[p]+1==dep[e[i].fi]){
				ll f=dfs(e[i].fi,min((ll)e[i].se,lim-sum));
				e[i].se-=f,e[i^1].se+=f;
				if((sum+=f)==lim)break;
			}
			if(!sum)dep[p]=0;
			return sum;
		};
		ll ret=0;
		while(bfs())ret+=dfs(s,infl);
		return ret;
	}
    inline V<V<pil>> gomoryhu(){
        V<int>f(e.size()),id(n),tmp(n);
        iota(ALL(id),0);
        shuffle(ALL(id),mt19937(time(0)));
        V<V<pil>>to(n);
        auto build=[&](auto &&self,int l,int r)->void{
            if(l==r)return;
            For(i,e.size())f[i]=e[i].se;
            ll w=dinic(id[l],id[r]);
            to[id[l]].eb(id[r],w),to[id[r]].eb(id[l],w);
            int L=l,R=r;
            FOR(i,l,r+1){
                if(dep[id[i]])tmp[L++]=id[i];
                else tmp[R--]=id[i];
            }
            assert(L==R+1);
            copy(tmp.begin()+l,tmp.begin()+r+1,id.begin()+l);
            For(i,e.size())e[i].se=f[i];
            self(self,l,R),self(self,L,r);
        };
        build(build,0,n-1);
        return to;
    }
};

inline pair<ll,V<int>> boundflow(int n,const V<array<int,4>> &e,int tp=0,int S=-1,int T=-1){
    // 0/1/2 can/max/min
    assert(tp>=0),assert(tp<=2);
    if(tp)assert(S>=0),assert(S<n),assert(T>=0),assert(T<n);
    else assert(!~S==!~T);
    if(S==T)S=T=-1,tp=0;
    V<ll>d(n);
    maxflow mf(n+2,n,n+1);
    for(const auto &[u,v,l,r]:e){
        assert(u>=0),assert(u<n),assert(v>=0),assert(v<n);
        assert(l>=0),assert(l<=r);
        d[u]-=l,d[v]+=l;
        mf.add_edge(u,v,r-l);
    }
    if(~S)mf.add_edge(T,S,infl);
    ll sum=0;
    For(i,n){
        if(d[i]>0)mf.add_edge(mf.S,i,d[i]),sum+=d[i];
        if(d[i]<0)mf.add_edge(i,mf.T,-d[i]);
    }
    ll f=mf.dinic();
    assert(f<=sum);
    if(f<sum)return {-1,{}};
    if(~T){
        Rep(i,n)if(d[i])mf.e.qb(),mf.e.qb(),mf.hd[i].qb();
        f=mf.e[mf.hd[S].back()].se;
        mf.e.qb(),mf.e.qb(),mf.hd[S].qb(),mf.hd[T].qb();
        mf.hd.qb(),mf.hd.qb(),mf.n=n;       
    }
    else f=0;
    if(tp&1)f+=mf.dinic(S,T);
    if(tp&2)ckmax(f-=mf.dinic(T,S),0ll); // 流量平衡时残量网络等于原网络，可能导致答案为负
    V<int>p(n),ret;
    ret.reserve(e.size());
    for(const auto &[u,v,l,r]:e)ret.pb(l+mf.e[mf.hd[v][p[v]++]].se),++p[u];
    return {f,ret};
}

struct mincost{
	V<array<int,3>>e;
	V<V<int>>hd;
	int n,S,T;
	inline void add_edge(int x,int y,int z,int w){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(z>=0);
		hd[x].pb(e.size()),e.pb({y,z,w}),hd[y].pb(e.size()),e.pb({x,0,-w});
	}
	inline mincost(int _n=0,int _S=-1,int _T=-1){
		V<V<int>>(n=_n).swap(hd);
		S=_S,T=_T;
	}
	inline mincost(const V<V<array<int,3>>> &to,int _S=-1,int _T=-1){
        V<V<int>>(n=to.size()).swap(hd);
		For(i,n)for(const array<int,3> &j:to[i])add_edge(i,j[0],j[1],j[2]);
		S=_S,T=_T;
	}
    typedef pair<ll,ll> pll;
	inline pll primal_dual(){
		assert(S!=-1),assert(T!=-1);
		V<ll>h;
        V<bool>vis(n);
        auto spfa=[&]{
            h.assign(n,infl);
            h[S]=0;
            queue<int>q;
            q.push(S);
            while(q.size()){
                int p=q.front();q.pop();
                vis[p]=false;
                for(int i:hd[p])if(e[i][1]&&ckmin(h[e[i][0]],h[p]+e[i][2])&&!vis[e[i][0]])q.push(e[i][0]),vis[e[i][0]]=true;
            }
        };
        spfa();
        V<ll>dis;
        V<pii>pre(n);
		auto dijkstra=[&](){
			V<ll>(n,infl).swap(dis);
            dis[S]=0;
			priority_queue<pli>q;
            q.emplace(0,S);
            V<bool>vis(n);
			while(q.size()){
				int p=q.top().se;q.pop();
                if(vis[p])continue;
                vis[p]=true;
				for(int i:hd[p])if(e[i][1]&&ckmin(dis[e[i][0]],dis[p]+e[i][2]+h[p]-h[e[i][0]])){
                    pre[e[i][0]]={p,i};
                    if(!vis[e[i][0]])q.emplace(-dis[e[i][0]],e[i][0]);
                }
			}
			return dis[T]!=infl;
		};
		ll ret1=0,ret2=0;
		while(dijkstra()){
            For(i,n)h[i]+=dis[i];
            ll f=infl;
            for(int i=T;i!=S;i=pre[i].fi)ckmin(f,(ll)e[pre[i].se][1]);
            for(int i=T;i!=S;i=pre[i].fi)e[pre[i].se][1]-=f,e[pre[i].se^1][1]+=f;
            ret1+=f,ret2+=f*h[T];
        }
		return {ret1,ret2};
	}
    inline pll dinic(){
        assert(S!=-1),assert(T!=-1);
        V<int>cur(n);
        V<ll>dis;
        V<V<int>>tmp=hd;
        V<bool>vis(n);
        auto spfa=[&](){
            dis.assign(n,infl);
            dis[S]=0;
            hd=tmp;
            queue<int>q;
            q.push(S);
            while(q.size()){
                int p=q.front();q.pop();
                vis[p]=false;
                for(int i:hd[p])if(e[i][1]&&ckmin(dis[e[i][0]],dis[p]+e[i][2])&&!vis[e[i][0]])q.push(e[i][0]),vis[e[i][0]]=true;
            }
            return dis[T]<infl;
        };
        ll ret1=0,ret2=0;
        auto dfs=[&](auto &&self,int p,ll f)->ll{
            if(p==T)return f;
            vis[p]=true;
            ll ret=0;
            while(hd[p].size()){
                int i=hd[p].back();
                if(!vis[e[i][0]]&&e[i][1]&&dis[e[i][0]]==dis[p]+e[i][2]){
                    ll d=self(self,e[i][0],min((ll)e[i][1],f-ret));
                    if(d){
                        ret+=d,ret2+=d*e[i][2];
                        e[i][1]-=d,e[i^1][1]+=d;
                        if(ret==f)break;
                    }
                }
                hd[p].qb();
            }
            vis[p]=false;
            return ret;
        };
        while(spfa()){
            ll d;
            while(d=dfs(dfs,S,infl))ret1+=d;
        }
        return {ret1,ret2};
    }
};
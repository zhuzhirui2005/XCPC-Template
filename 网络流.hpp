struct maxflow{
	V<pii>e;
	V<V<int>>hd;
	int n,S,T;
	inline void add_edge(int x,int y,int z){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(z>=0);
		hd[x].pb(e.size()),e.eb(y,z),hd[y].pb(e.size()),e.eb(x,0);
	}
	inline maxflow(int _n=0,int _S=-1,int _T=-1){
		V<V<int>>(n=_n).swap(hd);
		S=_S,T=_T;
	}
	inline maxflow(const V<V<pii>> &to,int _S=-1,int _T=-1){
        V<V<int>>(n=to.size()).swap(hd);
		For(i,n)for(const pii &j:to[i])add_edge(i,j.fi,j.se);
		S=_S,T=_T;
	}
	inline ll dinic(){
		assert(S!=-1),assert(T!=-1);
		V<int>dep;
		auto bfs=[&](){
			V<int>(n).swap(dep);
			dep[S]=1;
			queue<int>q;
			q.push(S);
			while(q.size()){
				int p=q.front();q.pop();
				for(int i:hd[p])if(e[i].se&&!dep[e[i].fi])dep[e[i].fi]=dep[p]+1,q.push(e[i].fi);
			}
			return dep[T];
		};
		function<ll(int,ll)>dfs=[&](int p,ll lim){
			if(p==T)return lim;
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
		while(bfs())ret+=dfs(S,infl);
		return ret;
	}
};

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
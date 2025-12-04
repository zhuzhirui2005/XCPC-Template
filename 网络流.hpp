struct maxflow{
    V<int>dep;
	V<pil>e;
	V<V<int>>hd;
	int n,S,T;
	inline void add_edge(int x,int y,ll z){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(z>=0);
		hd[x].pb(e.size()),e.eb(y,z),hd[y].pb(e.size()),e.eb(x,0);
	}
	inline maxflow(int n=0,int S=-1,int T=-1):n(n),S(S),T(T),hd(n){}
    template<class U>inline maxflow(const V<V<pair<int,U>>> &to,int S=-1,int T=-1):n(to.size()),S(S),T(T),hd(to.size()){For(i,n)for(const auto &[j,k]:to[i])add_edge(i,j,k);}
	inline ll dinic(int s=-1,int t=-1){
        if(!~s)s=S;
        if(!~t)t=T;
		assert(~s),assert(~t);
        queue<int>q;
		auto bfs=[&](){
			dep.assign(n,0);
			dep[s]=1;
			q.push(s);
			while(q.size()){
				int p=q.front();q.pop();
				for(int i:hd[p]){
                    const auto &[j,k]=e[i];
                    if(k&&!dep[j])dep[j]=dep[p]+1,q.push(j);
                }
			}
			return dep[t];
		};
        V<int>cur;
		auto dfs=[&](auto &&self,int p,ll lim){
			if(p==t)return lim;
			ll sum=0;
			for(int &idx=cur[p];idx<hd[p].size();++idx){
                int i=hd[p][idx];
                auto &[j,k]=e[i];
                if(k&&dep[p]+1==dep[j]){
                    ll f=self(self,j,min((ll)k,lim-sum));
                    k-=f,e[i^1].se+=f;
                    if((sum+=f)==lim)break;
                }
            }
			if(!sum)dep[p]=0;
			return sum;
		};
		ll ret=0;
		while(bfs())cur.assign(n,0),ret+=dfs(dfs,s,infl);
		return ret;
	}
    inline V<V<pil>> gomoryhu(){
        V<ll>f(e.size());
        V<int>id(n),tmp(n);
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
    for(const auto &[u,v,l,r]:e)ret.pb(l+mf.e[mf.hd[v][u==v?++p[v]:p[v]++]].se),++p[u];
    return {f,ret};
}

struct mincost{
    V<ll>dis;
	V<array<int,3>>e;
	V<V<int>>hd;
	int n,S,T;
	inline void add_edge(int x,int y,int z,int w){
		assert(0<=x),assert(x<n),assert(0<=y),assert(y<n),assert(z>=0);
		hd[x].pb(e.size()),e.pb({y,z,w}),hd[y].pb(e.size()),e.pb({x,0,-w});
	}
	inline mincost(int n=0,int S=-1,int T=-1):n(n),S(S),T(T),hd(n){}
	inline mincost(const V<V<array<int,3>>> &to,int S=-1,int T=-1):n(to.size()),S(S),T(T),hd(to.size()){For(i,n)for(const array<int,3> &j:to[i])add_edge(i,j[0],j[1],j[2]);}
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
                for(int i:hd[p]){
                    const auto &[j,k,l]=e[i];
                    if(k&&ckmin(h[j],h[p]+l)&&!vis[j])q.push(j),vis[j]=true;
                }
            }
        };
        spfa();
        V<pii>pre(n);
        priority_queue<pli>q;
		auto dijkstra=[&](){
			dis.assign(n,infl);
            dis[S]=0;
            q.emplace(0,S);
            vis.assign(n,false);
			while(q.size()){
				int p=q.top().se;q.pop();
                if(vis[p])continue;
                vis[p]=true;
				for(int i:hd[p]){
                    const auto &[j,k,l]=e[i];
                    if(k&&ckmin(dis[j],dis[p]+l+h[p]-h[j])){
                        pre[j]={p,i};
                        if(!vis[j])q.emplace(-dis[j],j);
                    }
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
        queue<int>q;
        V<bool>vis(n);
        auto spfa=[&](){
            dis.assign(n,infl);
            dis[S]=0;
            q.push(S);
            while(q.size()){
                int p=q.front();q.pop();
                vis[p]=false;
                for(int i:hd[p]){
                    const auto &[j,k,l]=e[i];
                    if(k&&ckmin(dis[j],dis[p]+l)&&!vis[j])q.push(j),vis[j]=true;
                }
            }
            return dis[T]<infl;
        };
        V<int>cur;
        ll ret1=0,ret2=0;
        auto dfs=[&](auto &&self,int p,ll f)->ll{
            if(p==T)return f;
            vis[p]=true;
            ll ret=0;
            for(int &idx=cur[p];idx<hd[p].size();++idx){
                int i=hd[p][idx];
                auto &[j,k,l]=e[i];
                if(!vis[j]&&k&&dis[j]==dis[p]+l){
                    ll d=self(self,j,min((ll)k,f-ret));
                    if(d){
                        ret+=d,ret2+=d*l;
                        k-=d,e[i^1][1]+=d;
                        if(ret==f)break;
                    }
                }
            }
            vis[p]=false;
            return ret;
        };
        while(spfa()){
            cur.assign(n,0);
            ll d;
            while(d=dfs(dfs,S,infl))ret1+=d;
        }
        return {ret1,ret2};
    }
};

struct maxmatch{
    V<int>dep,mch;
    int m,n;
    V<V<int>>to;
    inline void add_edge(int x,int y){
        assert(0<=x),assert(x<n),assert(0<=y),assert(y<m);
        to[x].pb(y);
    }
    inline maxmatch(int n=0,int m=0):n(n),m(m),mch(n+m,-1),to(n){}
    inline maxmatch(const V<V<int>> &to):n(to.size()),m(0),to(to){
        For(i,n)for(int j:to[i])assert(j>=0),ckmax(m,j);
        ++m;
        mch.assign(n+m,-1);
    }
    inline int hopcroft_karp(){
        dep.resize(n);
        auto bfs=[&](){
            queue<int>q;
            For(i,n){
                if(~mch[i])dep[i]=0;
                else dep[i]=1,q.push(i);
            }
            int mn=inf;
            while(q.size()){
                int p=q.front();q.pop();
                if(dep[p]>=mn)break;
                for(int i:to[p]){
                    int j=mch[n+i];
                    if(!~j)mn=dep[p]+1;
                    else if(!dep[j])dep[j]=dep[p]+1,q.push(j);
                }
            }
            return mn<inf;
        };
        V<int>cur;
        auto dfs=[&](auto &&self,int p)->bool{
            for(int &idx=cur[p];idx<to[p].size();++idx){
                int i=to[p][idx],&j=mch[n+i];
                if(!~j||(dep[j]==dep[p]+1&&self(self,j))){
                    mch[p]=n+i,j=p;
                    return true;
                }
            }
            dep[p]=0;
            return false;
        };
        int cnt=0;
        while(bfs()){
            cur.assign(n,0);
            For(i,n)if(!~mch[i]&&dfs(dfs,i))++cnt;
        }
        return cnt;
    }
};
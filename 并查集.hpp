struct dsu{
    V<int>fa;
    inline void resize(int n){V<int>(n,-1).swap(fa);}
    inline dsu(int n=0){resize(n);}
    int find(int k){return fa[k]<0?k:fa[k]=find(fa[k]);}
    inline bool merge(int x,int y){
		x=find(x),y=find(y);
		if(x!=y)fa[x]+=fa[y],fa[y]=x;
		return x!=y;
	}
	inline bool same(int x,int y){return find(x)==find(y);}
	inline int size(int k){return -fa[find(k)];}
};
inline pair<V<V<int>>,V<int>> kruskal_tree(int n,V<array<int,3>> &e){
    int cnt=n;
    dsu d(n+n-1);
    V<V<int>>to(n+n-1);
    V<int>val(n+n-1);
    sort(ALL(e),[&](const auto &x,const auto &y){return x[2]<y[2];});
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

struct range_dsu{
	V<V<int>>fa;
	int lg,n;
    inline void resize(int _n){
		V<V<int>>(lg=((n=_n)?__lg(n):-1)+1).swap(fa);
		For(i,lg)fa[i].resize(n-(1<<i)+1,-1);
	}
    inline range_dsu(int _n=0){resize(_n);}
    int find(int d,int k){return fa[d][k]<0?k:fa[d][k]=find(d,fa[d][k]);}
    inline void merge(int d,int x,int y){
		x=find(d,x),y=find(d,y);
		if(x>y)swap(x,y);
		if(x!=y)fa[d][x]+=fa[d][y],fa[d][y]=x;
	}
    inline void merge(int x1,int x2,int y1,int y2){
    	assert(x2-x1==y2-y1);
		Rep(i,lg)if(x1+(1<<i)-1<=x2){
			merge(i,x1,y1);
			x1+=1<<i,y1+=1<<i;
		}
	}
	inline void init(){
		REP(i,1,lg)For(j,n-(1<<i)+1){
			int k=find(i,j);
			merge(i-1,j,k),merge(i-1,j+(1<<i-1),k+(1<<i-1));
		}
	}
	int find(int k){return fa[0][k]<0?k:fa[0][k]=find(fa[0][k]);}
	inline bool same(int x,int y){return find(x)==find(y);}
	inline int size(int k){return -fa[0][find(k)];}
};
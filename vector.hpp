template<class T>
inline int area_cnt(const V<V<T>>& v,const T& goal){
    int n=v.size(),m=v[0].size(),ret=0;
    V<V<bool>>vis(n,V<bool>(m));
    function<void(int,int)>dfs=[&](int x,int y){
        vis[x][y]=true;
        For(i,4){
            int mx=x+dir[i][0],my=y+dir[i][1];
            if(mx>=0&&mx<n&&my>=0&&my<m&&!vis[mx][my]&&v[mx][my]==goal)dfs(mx,my);
        }
    };
    For(i,n)For(j,m)if(!vis[i][j]&&v[i][j]==goal)dfs(i,j),++ret;
    return ret;
}

inline array<V<V<int>>,4> find_dir(const V<V<int>>& grid){
	// Used "array<V<V<int>>,4>" in order to use "auto [u,d,l,r]=find(grid);"
	int n=grid.size(),m=grid[0].size();
	V<V<int>>u(n,V<int>(m)),d(n,V<int>(m)),l(n,V<int>(m)),r(n,V<int>(m));
	For(i,n)
	    For(j,m)
	        if(grid[i][j]){
	            u[i][j]=(i?u[i-1][j]:0)+1;
	            l[i][j]=(j?l[i][j-1]:0)+1;
	        }
	Rep(i,n)
	    Rep(j,m)
	        if(grid[i][j]){
	            d[i][j]=(i+1<n?d[i+1][j]:0)+1;
	            r[i][j]=(j+1<m?r[i][j+1]:0)+1;
	        }
	// Up   , Down , Left, Right
	// North, South, West, East
	return {u,d,l,r};
}

template<class T>
inline V<V<T>> rot(const V<V<T>>& v){
    V<V<T>>ret(v[0].size(),V<T>(v.size()));
    For(i,v.size())
        For(j,v[0].size())
            ret[j][v.size()-i-1]=v[i][j];
    return ret;
}
template<class T>
inline V<V<T>> rot(const V<V<T>>&v,int k){
    assert(k>=0);
    k&=3;
    auto ret=v;
    while(k--)ret=rot(ret);
    return ret;
}

struct ListNode {
     int val;
     ListNode *next;
     ListNode() : val(0), next(nullptr) {}
     ListNode(int x) : val(x), next(nullptr) {}
     ListNode(int x, ListNode *next) : val(x), next(next) {}
};
inline V<int> ltov(ListNode *hd){
	V<int>ret;
	while(hd)ret.pb(hd->val),hd=hd->next;
	return ret;
}
inline ListNode* vtol(const V<int>& v){
	if(v.empty())return NULL;
	ListNode *hd=new ListNode(v[0]),*p=hd;
	FOR(i,1,v.size())p->next=new ListNode(v[i]),p=p->next;
	return hd;
}

template<class T1,class T2>
inline V<T1> pre1d(const V<T2> &v,const function<T1(T1,T2)> &f=[](const T1 &x,const T2 &y){return x+y;}){
    V<T1>ret(v.size());
    For(i,v.size())ret[i]=f((i?ret[i-1]:T1()),v[i]);
    return ret;
}
template<class T1,class T2>
inline void add1d(V<T1> &v,int l,int r,const T2 &k,const function<T1(T1,T2)> &f=[](const T1 &x,const T2 &y){return x+y;}){
    assert(0<=l),assert(l<=r),assert(r<v.size());
    v[l]=f(v[l],k);
    if(r+1<v.size())v[r+1]=f(v[r+1],-k);
}
template<class T>
inline T qry1d(V<T> &v,int l,int r,const function<T(T,T)> &f=[](const T &x,const T &y){return x+y;}){
    assert(0<=l),assert(l<=r),assert(r<v.size());
    return f(v[r],-(l?v[l-1]:T()));
}

template<class T1,class T2>
inline V<V<T1>> pre2d(const V<V<T2>> &v){
    V<V<T1>>ret(v.size(),V<T1>(v[0].size()));
    For(i,v.size())For(j,v[0].size())ret[i][j]=v[i][j]+(i?ret[i-1][j]:T1())+(j?ret[i][j-1]:T1())-(i&&j?ret[i-1][j-1]:T1());
    return ret;
}
template<class T1,class T2>
inline void add2d(V<V<T1>> &v,int l1,int l2,int r1,int r2,const T2 &k){
    assert(0<=l1),assert(l1<=r1),assert(r1<v.size());
    assert(0<=l2),assert(l2<=r2),assert(r2<v[0].size());
    v[l1][l2]+=k;
    if(r1+1<v.size())v[r1+1][l2]-=k;
    if(r2+1<v[0].size())v[l1][r2+1]-=k;
    if(r1+1<v.size()&&r2+1<v[0].size())v[r1+1][r2+1]+=k;
}
template<class T>
inline T qry2d(V<V<T>> &v,int l1,int l2,int r1,int r2){
    assert(0<=l1),assert(l1<=r1),assert(r1<v.size());
    assert(0<=l2),assert(l2<=r2),assert(r2<v[0].size());
    return v[r1][r2]-(l1?v[l1-1][r2]:T())-(l2?v[r1][l2-1]:T())+(l1&&l2?v[l1-1][l2-1]:T());
}

inline ll contor(const V<int> &v){
	int d=*min_element(ALL(v)),n=v.size();
	V<bool>vis(n);
	for(int i:v)vis[i-d]=true;
	if(any_of(ALL(vis),[](bool b){return !b;}))return -1;
	V<ll>fac(n);
	fac[0]=1;
	BIT3<int>t(n);
	FOR(i,1,n+1){
		if(i<n)fac[i]=fac[i-1]*i;
		++t.c[i];
		if(i+(i&-i)<=n)t.c[i+(i&-i)]+=t.c[i];
	}
	ll ret=0;
	For(i,n){
		t.add(v[i]-d,-1);
		ret+=fac[n-i-1]*t.query(v[i]-d);
	}
	return ret;
}

inline V<int> inv_contor(int n,ll k){
	V<ll>fac(n+1);
	fac[0]=1;
	FOR(i,1,n+1)fac[i]=fac[i-1]*i;
	if(k>=fac[n])return {-1};
	V<int>ret(n);
	V<bool>vis(n);
	For(i,n){
		int dgt=k/fac[n-i-1]+1,j=-1;
		k%=fac[n-i-1];
		do dgt-=!vis[++j];while(dgt);
		ret[i]=j,vis[j]=true;
	}
	return ret;
}
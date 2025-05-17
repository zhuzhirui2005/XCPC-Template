template<class T>
inline V<V<T>> rot(const V<V<T>>& v){
    V<V<T>>ret(v[0].size(),V<T>(v.size()));
    For(i,v.size())
        For(j,v[0].size())
            ret[j][v.size()-i-1]=v[i][j];
    return ret;
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
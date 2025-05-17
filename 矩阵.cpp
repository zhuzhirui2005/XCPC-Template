template<class T>
struct matrix{
	int n,m;
	V<V<T>>a;
	inline matrix(int _n=0,int _m=0,T v=T()):n(_n),m(_m){
		V<V<T>>(n,V<T>(m,v)).swap(a);
	};
	inline V<T> &operator[](int idx){return a[idx];}
	inline const V<T> &operator[](int idx)const{return a[idx];}
	inline matrix operator*(const matrix &rhs){
		assert(m==rhs.n);
		matrix ret(n,rhs.m);
		For(i,n)For(j,rhs.m)For(k,m)ret[i][j]+=a[i][k]*rhs[k][j];
		return ret;
	}
	inline matrix trans(){
		matrix ret(m,n);
		For(i,n)For(j,m)ret[j][i]=a[i][j];
		return ret;
	}
	inline bool gauss(){
		assert(n<=m);
		int nw=0;
		For(i,n){
			if(!a[nw][i])FOR(j,nw+1,n)if(a[j][i]){swap(a[nw],a[j]);break;}
			if(a[nw][i]){
				For(j,n)if(nw!=j){
					T coef=a[j][i]/a[nw][i];
					FOR(k,i,m)a[j][k]-=coef*a[nw][k];
				}
				++nw;
			}
		}
		return nw==n;
	}
	inline matrix unit(){
		assert(n==m);
		matrix ret(n,n);
		For(i,n)ret[i][i]=1;
		return ret;
	}
	inline matrix pow(ull k){
		matrix base=*this,ret=unit(); 
		for(;k;k>>=1,base=base*base)if(k&1)ret=ret*base;
		return ret;
	}
	inline matrix mul_pow(const matrix &rhs,ull k){
		matrix base=rhs,ret=*this; 
		for(;k;k>>=1,base=base*base)if(k&1)ret=ret*base;
		return ret;
	}
};

template<class T>
struct dis_matrix{
	int n,m;
	V<V<T>>a;
	inline dis_matrix(int _n=0,int _m=0,T v=T()):n(_n),m(_m){
		assert((is_same<T,int>::value)||(is_same<T,ll>::value)||(is_same<T,ull>::value));
		V<V<T>>(n,V<T>(m)).swap(a);
	};
	inline V<T> &operator[](int idx){return a[idx];}
	inline const V<T> &operator[](int idx)const{return a[idx];}
	inline dis_matrix operator*(const dis_matrix &rhs){
		assert(m==rhs.n);
		dis_matrix ret(n,rhs.m,is_same<T,int>::value?inf:infl);
		For(i,n)For(j,rhs.m)For(k,m)ckmin(ret[i][j],a[i][k]+rhs[k][j]);
		return ret;
	}
	inline dis_matrix pow(ull k){
		dis_matrix base=*this,ret(n,n,is_same<T,int>::value?inf:infl);
		for(;k;k>>=1,base=base*base)if(k&1)ret=ret*base;
		return ret;
	}
};
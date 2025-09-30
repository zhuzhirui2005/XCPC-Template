template<class T,T(*merge)(T,T)>
struct ST{
	V<V<T>>st;
	inline ST(){}
	inline ST(const V<T> &a){
		int n=a.size(),B=__lg(n);
		V<V<T>>(B+1).swap(st);
		st[0]=a;
		FOR(i,1,B+1){
			st[i].resize(n-(1<<i)+1);
			For(j,n-(1<<i)+1)st[i][j]=merge(st[i-1][j],st[i-1][j+(1<<i-1)]);
		}
	}
	inline ST(const V<T> &a,const V<int> &pos){
		assert(a.size()==pos.size());
		int n=a.size(),B=__lg(n);
		V<V<T>>(B+1).swap(st);
		For(i,B+1){
			st[i].resize(n-(1<<i)+1);
			if(i)For(j,n-(1<<i)+1)st[i][j]=merge(st[i-1][j],st[i-1][j+(1<<i-1)]);
			else For(i,n)st[0][pos[i]]=a[i];
		}
	}
	inline T query(int l,int r){
		int n=st[0].size();
		assert(0<=l),assert(l<=r),assert(r<n);
		int k=__lg(r-l+1);
		return merge(st[k][l],st[k][r-(1<<k)+1]);
	}
};
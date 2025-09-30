template<class T,int n>
struct LB{
	V<T>d;
	int cnt,failed;
	inline void clear(){cnt=failed=0,V<T>(n).swap(d);}
	inline LB(){
		assert(n>0);
		assert(n<=(is_same<T,int>::value?31:is_same<T,unsigned>::value?32:is_same<T,ll>::value?63:is_same<T,ull>::value?64:-1));
		clear();
	}
	inline bool insert(T k){
		Rep(i,n)if(k>>i&1){
			if(!d[i]){
				++cnt,d[i]=k;
				return true;
			}
			else if(!(k^=d[i]))break;
		}
		++failed;
		return false;
	}
	inline bool can(T k){
		Rep(i,n)if(k>>i&1){
			if(!d[i])return false;
			else if(!(k^=d[i]))break;
		}
		return true;
	}
	inline T mx(T k=0){
		Rep(i,n)ckmax(k,k^d[i]);
		return k;
	}
	inline LB operator+(const LB &rhs){
		LB ret=rhs;
		ret.failed+=failed;
		For(i,n)if(d[i])ret.insert(d[i]);
		return ret;
	}
	inline LB &operator+=(const LB &rhs){
		failed+=rhs.failed;
		For(i,n)if(rhs.d[i])insert(rhs.d[i]);
		return *this;
	}
	// assumed that empty set isn't allowed
	inline T count(){return (T(1)<<cnt)-!failed;}
	// 0-indexed
	inline T rk(T k){
		T pw2=1,ret=0;
		For(i,n)if(d[i]){
			if(k>>i&1)ret|=pw2;
			pw2<<=1;
		}
		return ret;
	}
	inline T at(T k){
		if(!failed)++k;
		FOR(i,1,n)Rep(j,i)if(d[i]>>j&1)d[i]^=d[j];
		T ret=0;
		For(i,n)if(d[i]){
			if(k&1)ret^=d[i];
			k>>=1;
		}
		return k?-1:ret;
	}
};

template<class T,int n>
struct LB_ts{ // timestamp
	V<T>d;
	V<int>t;
	inline void clear(){V<T>(n).swap(d),V<int>(n).swap(t);}
	inline LB_ts(){
		assert(n>0);
		assert(n<=(is_same<T,int>::value?31:is_same<T,unsigned>::value?32:is_same<T,ll>::value?63:is_same<T,ull>::value?64:-1));
		clear();
	}
	inline bool insert(T k,int tm){
		Rep(i,n)if(k>>i&1){
			if(!d[i]){
				d[i]=k,t[i]=tm;
				return true;
			}
			else if(tm>t[i])swap(d[i],k),swap(t[i],tm);
			if(!(k^=d[i]))break;
		}
		return false;
	}
	inline bool can(T k,int tm=0){
		Rep(i,n)if(k>>i&1){
			if(!d[i]||t[i]<tm)return false;
			else if(!(k^=d[i]))break;
		}
		return true;
	}
	inline T mx(T k=0,int tm=0){
		Rep(i,n)if(t[i]>=tm)ckmax(k,k^d[i]);
		return k;
	}
	inline LB_ts operator+(const LB_ts &rhs){
		LB_ts ret=rhs;
		For(i,n)if(d[i])ret.insert(d[i],t[i]);
		return ret;
	}
	inline LB_ts &operator+=(const LB_ts &rhs){
		For(i,n)if(rhs.d[i])insert(rhs.d[i],rhs.t[i]);
		return *this;
	}
};
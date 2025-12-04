template<class T,class U=less<T>>
struct ST{
    const U cmp=U();
	int n;
	V<V<T>>st;
	inline ST(){}
	inline ST(const V<T> &a):n(a.size()){
        if(n){
            int B=__lg(n);
            V<V<T>>(B+1).swap(st);
            st[0]=a;
            FOR(i,1,B+1){
                st[i].resize(n-(1<<i)+1);
                For(j,n-(1<<i)+1)st[i][j]=min(st[i-1][j],st[i-1][j+(1<<i-1)],cmp);
            }
        }
	}
    inline ST(const V<T> &a,const V<int> &pos):n(a.size()){
		assert(n==pos.size());
		if(n){
			int B=__lg(n);
			V<V<T>>(B+1).swap(st);
			For(i,B+1){
				st[i].resize(n-(1<<i)+1);
				if(i)For(j,n-(1<<i)+1)st[i][j]=min(st[i-1][j],st[i-1][j+(1<<i-1)],cmp);
				else For(i,n)st[0][pos[i]]=a[i];
			}
		}
	}
	inline T query(int l,int r){
		assert(0<=l),assert(l<=r),assert(r<n);
		int k=__lg(r-l+1);
		return min(st[k][l],st[k][r-(1<<k)+1],cmp);
	}
};

template<class T,class U=less<T>>
struct RMQ{
    V<T>a,pre,suf;
    const U cmp=U();
    int n;
	V<V<T>>st;
    V<ull>stk;
    inline RMQ(){}
    inline RMQ(const V<T> &a):a(a),n(a.size()),pre(a),suf(a),stk(a.size()){
        if(n){
            const int m=n+63>>6,B=__lg(m)+1;
            st.resize(B);
            For(i,B){
                st[i].resize(m-(1<<i)+1);
                if(i)For(j,m-(1<<i)+1)st[i][j]=min(st[i-1][j],st[i-1][j+(1<<i-1)],cmp);
                else For(j,m){
                    st[0][j]=a[j<<6];
                    FOR(k,1,min(64,n-(j<<6)))st[0][j]=min(st[0][j],a[j<<6|k],cmp);
                }
            }
            FOR(i,1,n)if(i&63)pre[i]=min(pre[i],pre[i-1],cmp);
            Rep(i,n-1)if((i&63)<63)suf[i]=min(suf[i],suf[i+1],cmp);
            For(i,m){
                ull S=0;
                FOR(j,i<<6,min(i+1<<6,n)){
                    while(S){
                        int k=__lg(S);
                        if(cmp(a[j],a[i<<6|k]))S^=1ull<<k;
                        else break;
                    }
                    stk[j]=(S|=1ull<<(j&63));
                }
            }
        }
    }
    inline T query(int l,int r){
        assert(0<=l),assert(l<=r),assert(r<n);
        int L=l>>6,R=r>>6;
        if(L<R){
            T ret=min(suf[l],pre[r],cmp);
            if(++L<=--R){
                int k=__lg(R-L+1);
                ret=min({ret,st[k][L],st[k][R-(1<<k)+1]},cmp);
            }
            return ret;
        }
        else return a[__builtin_ctzll(stk[r]>>(l&63))+l];
    }
};
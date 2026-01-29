template<class T>
struct matrix{
    int n,m;
    V<V<T>>a;
    inline matrix(int n=0,int m=0,T v=T()):n(n),m(m),a(n,V<T>(m,v)){}
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
            if(!a[nw][i])FOR(j,nw+1,n)if(a[j][i]!=0){swap(a[nw],a[j]);break;}
            if(a[nw][i]){
                T inv=1/a[nw][i];
                For(j,n)if(nw!=j){
                    T coef=a[j][i]*inv;
                    FOR(k,i,m)a[j][k]-=coef*a[nw][k];
                }
                ++nw;
            }
        }
        return nw==n;
    }
    inline T det(){
        assert(n==m);
        T ret=1;
        For(i,n){
            if(!a[i][i])FOR(j,i+1,n)if(a[j][i]!=0){
                ret=-ret;
                swap(a[i],a[j]);
                break;
            }
            if(!a[i][i])return 0;
            T inv=1/a[i][i];
            ret*=a[i][i];
            FOR(j,i+1,n){
                T coef=a[j][i]*inv;
                FOR(k,i,m)a[j][k]-=a[i][k]*coef;
            }
        }
        return ret;
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
    inline matrix pow_mul(const matrix &rhs,ull k)const{
        matrix base=*this,ret=rhs; 
        for(;k;k>>=1,base=base*base)if(k&1)ret=base*ret;
        return ret;
    }
};

template<class T>
struct dis_matrix{
    static_assert((is_same_v<T,int>)||(is_same_v<T,ll>)||(is_same_v<T,ull>));
    int n,m;
    V<V<T>>a;
    inline dis_matrix(int n=0,int m=0,T v=T()):n(n),m(m),a(n,V<T>(m,v)){}
    inline V<T> &operator[](int idx){return a[idx];}
    inline const V<T> &operator[](int idx)const{return a[idx];}
    inline dis_matrix operator*(const dis_matrix &rhs){
        assert(m==rhs.n);
        dis_matrix ret(n,rhs.m,is_same_v<T,int>?inf:infl);
        For(i,n)For(j,rhs.m)For(k,m)ckmin(ret[i][j],a[i][k]+rhs[k][j]);
        return ret;
    }
    inline dis_matrix pow(ull k){
        dis_matrix base=*this,ret(n,n,is_same_v<T,int>?inf:infl);
        for(;k;k>>=1,base=base*base)if(k&1)ret=ret*base;
        return ret;
    }
    inline dis_matrix mul_pow(const dis_matrix &rhs,ull k){
        dis_matrix base=rhs,ret=*this; 
        for(;k;k>>=1,base=base*base)if(k&1)ret=ret*base;
        return ret;
    }
};
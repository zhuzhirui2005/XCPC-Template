template<class T,int n>
struct LB{
    static_assert(n>0);
    static_assert(n<=(is_same_v<T,int>?31:is_same_v<T,unsigned>?32:is_same_v<T,ll>?63:is_same_v<T,ull>?64:-1));
    int cnt,failed;
    array<T,n>d;
    inline LB():cnt(0),failed(0),d{}{}
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
        Rep(i,n)if(!(k>>i&1))k^=d[i];
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
    inline auto count(){
        // 0 means empty if failed is false, else full
        return (cnt?make_unsigned_t<T>(2)<<cnt-1:1)-!failed;
    }
    // 0-indexed
    inline T rk(T k){
        T pw2=1,ret=0;
        For(i,n)if(d[i]){
            if(k>>i&1)ret|=pw2;
            pw2<<=1;
        }
        return ret-!failed;
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
    static_assert(n>0);
    static_assert(n<=(is_same_v<T,int>?31:is_same_v<T,unsigned>?32:is_same_v<T,ll>?63:is_same_v<T,ull>?64:-1));
    array<T,n>d;
    array<int,n>t;
    inline LB_ts():d{},t{}{}
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
        Rep(i,n)if(t[i]>=tm&&!(k>>i&1))k^=d[i];
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
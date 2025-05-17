template<int p>
struct modint{
    int val;
    inline modint(int v=0):val(v){}
    inline modint& operator=(int v){val=v;return *this;}
    inline modint& operator+=(const modint&k){val=val+k.val>=p?val+k.val-p:val+k.val;return *this;}
    inline modint& operator-=(const modint&k){val=val-k.val<0?val-k.val+p:val-k.val;return *this;}
    inline modint& operator*=(const modint&k){val=int(1ll*val*k.val%p);return *this;}
    inline modint& operator^=(int k){modint r(1),b=*this;for(;k;k>>=1,b*=b)if(k&1)r*=b;val=r.val;return *this;}
    inline modint& operator/=(modint k){return *this*=(k^=p-2);}
    inline modint& operator+=(int k){val=val+k>=p?val+k-p:val+k;return *this;}
    inline modint& operator-=(int k){val=val<k?val-k+p:val-k;return *this;}
    inline modint& operator*=(int k){val=int(1ll*val*k%p);return *this;}
    inline modint& operator/=(int k){return *this*=((modint(k))^=p-2);}
    template<class T>friend modint operator+(modint a,T b){return a+=b;}
    template<class T>friend modint operator-(modint a,T b){return a-=b;}
    template<class T>friend modint operator*(modint a,T b){return a*=b;}
    template<class T>friend modint operator/(modint a,T b){return a/=b;}
    friend modint operator^(modint a,int b){return a^=b;}
    friend bool operator==(modint a,int b){return a.val==b;}
    friend bool operator!=(modint a,int b){return a.val!=b;}
    inline bool operator!()const{return !val;}
    inline modint operator-()const{return val?modint(p-val):modint(0);}
    inline modint operator++(int){modint t=*this;*this+=1;return t;}
    inline modint& operator++(){return *this+=1;}
    inline modint operator--(int){modint t=*this;*this-=1;return t;}
    inline modint& operator--(){return *this-=1;}
};
using mi=modint<mod>;
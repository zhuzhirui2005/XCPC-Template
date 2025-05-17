struct fraction{
    ll p,q;
    inline void simplify(){ll g=gcd(p<0?-p:p,q);p/=g;q/=g;}
    inline explicit fraction(ll _p=0):p(_p),q(1){}
    inline fraction(ll _p,ll _q):p(_p),q(_q){assert(q);if(q<0)p=-p,q=-q;simplify();}
    inline explicit fraction(const string&s){size_t pos=s.find('.');q=1;if(pos==string::npos)p=stoll(s);else{if(pos+1<s.size()){for(int i=0;i<s.size()-1-pos;i++)q*=10;p=(pos?stoll(s.substr(0,pos))*q:0)+stoll(s.substr(pos+1));}else p=stoll(s.substr(0,pos));simplify();}}
    inline explicit fraction(const V<char>&s):fraction(string(s.begin(),s.end())){}
    inline fraction& operator=(const fraction&r){p=r.p;q=r.q;return*this;}
    inline fraction& operator=(ll r){p=r;q=1;return*this;}
    inline fraction operator+(const fraction&r)const{if(q==r.q)return{p+r.p,q};ll g=gcd(q,r.q),m=q/g;return{p*(r.q/g)+r.p*m,m*r.q};}
    inline fraction operator+(ll r)const{return{p+r*q,q};}
    inline fraction add(const fraction&r)const{return{p*r.q+r.p*q,q*r.q};}
    inline fraction operator-(const fraction&r)const{if(q==r.q)return{p-r.p,q};ll g=gcd(q,r.q),m=q/g;return{p*(r.q/g)-r.p*m,m*r.q};}
    inline fraction operator-(ll r)const{return{p-r*q,q};}
    inline fraction sub(const fraction&r)const{return{p*r.q-r.p*q,q*r.q};}
    inline fraction operator*(const fraction&r)const{fraction t;ll g1=gcd(p,r.q),g2=gcd(r.p,q);t.p=(p/g1)*(r.p/g2);t.q=(q/g2)*(r.q/g1);return t;}
    inline fraction operator*(ll r)const{fraction t=*this;ll g=gcd(r,q);t.p*=r/g;t.q/=g;return t;}
    inline fraction mul(const fraction&r)const{return{p*r.p,q*r.q};}
    inline fraction operator/(const fraction&r)const{assert(r.p);fraction t;ll g1=gcd(p,r.p),g2=gcd(r.q,q);t.p=(p/g1)*(r.q/g2);t.q=(q/g2)*(r.p/g1);return t;}
    inline fraction operator/(ll r)const{assert(r);fraction t=*this;ll g=gcd(p,r);t.p/=g;t.q*=r/g;return t;}
    inline fraction div(const fraction&r)const{assert(r.p);return{p*r.q,q*r.p};}
    inline bool operator==(const fraction&r)const{return p==r.p&&q==r.q;}
    inline bool operator==(ll r)const{return p==r&&q==1;}
    inline bool eq(const fraction&r)const{return p==r.p&&q==r.q;}
    inline bool operator<(const fraction&r)const{return p*r.q<r.p*q;}
    inline bool operator<(ll r)const{ll g=gcd(p,r);return p/g<q*(r/g);}
    inline bool lt(const fraction&r)const{return p*r.q<r.p*q;}
    inline bool operator>(const fraction&r)const{return p*r.q>r.p*q;}
    inline bool operator>(ll r)const{ll g=gcd(p,r);return p/g>q*(r/g);}
    inline bool gt(const fraction&r)const{return p*r.q>r.p*q;}
    inline bool operator<=(const fraction&r)const{return p*r.q<=r.p*q;}
    inline bool operator<=(ll r)const{ll g=gcd(p,r);return p/g<=q*(r/g);}
    inline bool le(const fraction&r)const{return p*r.q<=r.p*q;}
    inline bool operator>=(const fraction&r)const{return p*r.q>=r.p*q;}
    inline bool operator>=(ll r)const{ll g=gcd(p,r);return p/g>=q*(r/g);}
    inline bool ge(const fraction&r)const{return p*r.q>=r.p*q;}
    inline string to_string()const{return ::to_string(p)+'/'+::to_string(q);}
};
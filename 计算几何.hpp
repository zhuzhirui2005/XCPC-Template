const double eps=1e-8,PI=acos(-1.);

struct vec2D{
    double x,y;
    inline vec2D(double x=0.,double y=0.):x(x),y(y){}
    inline vec2D operator+(const vec2D &rhs){return {x+rhs.x,y+rhs.y};}
    inline vec2D operator-(const vec2D &rhs)const{return {x-rhs.x,y-rhs.y};}
    inline vec2D operator*(double k){return {x*k,y*k};}
    inline double cross(const vec2D &rhs)const{return x*rhs.y-y*rhs.x;}
    inline double dot(const vec2D &rhs){return x*rhs.x+y*rhs.y;}
    inline double operator*(const vec2D &rhs){return dot(rhs);}
    inline bool operator==(const vec2D &rhs){return max(abs(x-rhs.x),abs(y-rhs.y))<eps;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2D &rhs){return dot(rhs)*dot(rhs)==dot(*this)*(rhs*rhs);}
    inline bool coln(const vec2D &rhs)const{return abs(cross(rhs))<eps;}
    inline int dir(const vec2D &rhs)const{return !coln(rhs)?cross(rhs)>=eps?1:-1:0;} // 1 if rhs is on ccw
    inline double dist(const vec2D &rhs){
        double dx=x-rhs.x,dy=y-rhs.y;
        return sqrt(dx*dx+dy*dy);
    }
    inline double norm(){return sqrt(x*x+y*y);}
    inline double projcoef(const vec2D &rhs){return dot(rhs)/dot(*this);}
    inline double projlen(const vec2D &rhs){return dot(rhs)/norm();}
    inline vec2D rot(double theta){ // ccw
        double c=cos(theta),s=sin(theta);
        return {x*c-y*s,y*c+x*s};
    }
    inline double ang(){ // ccw from positive x axis, same as vec2D(1,0).ang(*this)
        return atan2(y,x);
    }
    inline double ang(const vec2D &rhs){ // ccw rot
        double ret=atan2(cross(rhs),dot(rhs));
        if(ret<0)ret+=2*PI;
        return ret;
    }
    inline int quad()const{
        if(x>0&&y>=0)return 1;
        if(x<=0&&y>0)return 2;
        if(x<0&&y<=0)return 3;
        if(x>=0&&y<0)return 4;
        return 0;
    }
};

struct vec2d{
    ll x,y;
    inline vec2d(ll x=0,ll y=0):x(x),y(y){}
    inline vec2d operator+(const vec2d &rhs){return {x+rhs.x,y+rhs.y};}
    inline vec2d operator-(const vec2d &rhs)const{return {x-rhs.x,y-rhs.y};}
    inline vec2d operator*(ll k){return {x*k,y*k};}
    inline ll cross(const vec2d &rhs)const{return x*rhs.y-y*rhs.x;}
    inline ll dot(const vec2d &rhs){return x*rhs.x+y*rhs.y;}
    inline ll operator*(const vec2d &rhs){return dot(rhs);}
    inline bool operator==(const vec2d &rhs){return x==rhs.x&&y==rhs.y;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2D &rhs){return dot(rhs)*dot(rhs)==dot(*this)*(rhs*rhs);}
    inline bool coln(const vec2d &rhs)const{return !cross(rhs);}
    inline int dir(const vec2d &rhs)const{return cross(rhs)?cross(rhs)>0?1:-1:0;} // 1 if rhs is on ccw
    inline double dist(const vec2d &rhs){
        ll dx=x-rhs.x,dy=y-rhs.y;
        return sqrt(dx*dx+dy*dy);
    }
    inline double norm(){return sqrt(x*x+y*y);}
    inline double projcoef(const vec2d &rhs){return 1.*dot(rhs)/dot(*this);}
    inline double projlen(const vec2d &rhs){return dot(rhs)/norm();}
    inline vec2D rot(double theta){ // ccw, note that the result used double
        double c=cos(theta),s=sin(theta);
        return {x*c-y*s,y*c+x*s};
    }
    inline double ang(){ // ccw from positive x axis, same as vec2d(1,0).ang(*this)
        return atan2(y,x);
    }
    inline double ang(const vec2d &rhs){ // ccw rot
        double ret=atan2(cross(rhs),dot(rhs));
        if(ret<0)ret+=2*PI;
        return ret;
    }
    inline int quad()const{
        if(x>0&&y>=0)return 1;
        if(x<=0&&y>0)return 2;
        if(x<0&&y<=0)return 3;
        if(x>=0&&y<0)return 4;
        return 0;
    }
};

template<class vec>
inline V<vec2D> intersect(const vec &p1,const vec &p2,const vec &p3,const vec &p4){ // (p1,p2) and (p3,p4) are segments
    vec u=p2-p1,v=p4-p3,w=p3-p1;
    auto isz=[&](const vec2D &k){return abs(k.x)<eps&&abs(k.y)<eps;};
    if(isz({u.x,u.y})&&isz({v.x,v.y})){
        if(isz(w))return {{p1.x,p1.y}};
        return {};
    }
    auto on=[&](const vec &p1,const vec &p2,const vec &p3){return min(p1.x,p3.x)-eps<p2.x&&p2.x<max(p1.x,p3.x)+eps&&min(p1.y,p3.y)-eps<p2.y&&p2.y<max(p1.y,p3.y)+eps;};
    if(isz({u.x,u.y})){
        if(v.coln(w)&&on(p3,p1,p4))return {{p1.x,p1.y}};
        return {};
    }
    if(isz({v.x,v.y})){
        if(u.coln(w)&&on(p1,p3,p2))return {{p3.x,p3.y}};
        return {};
    }
    if(u.coln(v)){
        if(!v.coln(w))return {};
        V<vec2D>ret;
        auto check=[&](const vec &p){
            for(const auto &i:ret)if(isz({p.x-i.x,p.y-i.y}))return false;
            return true;
        };
        if(on(p3,p1,p4))ret.eb(p1.x,p1.y);
        if(on(p3,p2,p4)&&check(p2))ret.eb(p2.x,p2.y);
        if(on(p1,p3,p2)&&check(p3))ret.eb(p3.x,p3.y);
        if(on(p1,p4,p2)&&check(p4))ret.eb(p4.x,p4.y);
        return ret;
    }
    else{
        double denom=u.cross(v),t1=w.cross(v)/denom,t2=w.cross(u)/denom;
        if(-eps<t1&&t1<1+eps&&-eps<t2&&t2<1+eps)return {{p1.x+t1*u.x,p1.y+t1*u.y}};
        return {};
    }
}

template<class vec>
inline bool angle_cmp(const vec &x,const vec &y){
    return x.quad()!=y.quad()?x.quad()<y.quad():x.dir(y)==1; // eps may break sort, use atan2 instead
}

template<class vec,bool coln=false> // coln is always considered as false when using vec2D
struct convex{
    int n;
    V<vec>p;
    inline int nxt(int k)const{return k+1==p.size()?0:k+1;}
    inline int pre(int k){return k?k-1:p.size()-1;} // p should not be empty
    inline void add(const vec &k){++n,p.pb(k);}
    inline auto cross(const vec &a,const vec &b,const vec &c)const{
        // (b-a).cross(c-a)
        return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
    }
    inline void init(){
        sort(ALL(p),[&](const vec &u,const vec &v){return u.x!=v.x?u.x<v.x:u.y<v.y;});
        p.erase(unique(ALL(p)),p.end());
        if((n=p.size())<3)return;
        V<vec>hi,lo;
        auto bad=[&](const V<vec> &v,const vec &k){
            if(v.size()<2)return false;
            auto c=cross(v[v.size()-2],v.back(),k);
            return is_same_v<vec,vec2D>?c<eps:coln?c<0:c<=0;
        };
        For(i,n){
            while(bad(lo,p[i]))lo.qb();
            lo.pb(p[i]);
        }
        Rep(i,n){
            while(bad(hi,p[i]))hi.qb();
            hi.pb(p[i]);
        }
        lo.qb(),hi.qb();
        n=lo.size()+hi.size(),p.swap(lo),p.insert(p.end(),ALL(hi));
    }
    inline convex():n(0){}
    inline convex(const V<vec> &p):p(p){init();}
    inline double area(){
        if(n<3)return 0.;
        double ret=0.;
        For(i,n)ret+=p[i].cross(p[nxt(i)]);
        return ret*.5;
    }
    inline double peri(){
        if(n<2)return 0.;
        double ret=0.;
        For(i,n)ret+=p[i].dist(p[nxt(i)]);
        return ret;
    }
    inline double diam(){
        if(n<2)return 0.;
        if(n==2)return p[0].dist(p[1]);
        int j=1;
        double ret=0.;
        For(i,n){
            int k=nxt(i);
            while(true){
                int l=nxt(j);
                if(cross(p[i],p[k],p[j])<cross(p[i],p[k],p[l]))j=l;
                else break;
            }
            ckmax(ret,p[j].dist(p[i]));
            if(j!=k)ckmax(ret,p[j].dist(p[k]));
        }
        return ret;
    }
    // (-1,-1) means inside
    inline pii reach(const vec &k)requires is_same_v<vec,vec2D>{ // always allow walking on edges
        if(n<3)return {0,n-1};
        auto check=[&](int i){return cross(p[i],p[nxt(i)],k)<eps;};
        int neg=-1,pos=-1;
        (check(0)?pos:neg)=0,(check(n-1)?pos:neg)=n-1;
        if(!~pos){
            int l=1,r=n-2;
            while(l<=r){
                int mid=l+r>>1;
                if(check(mid)){pos=mid;break;}
                if(cross(p[mid],p[0],k)<eps)l=mid+1;
                else r=mid-1;
            }
            if(!~pos)return {-1,-1};
        }
        if(!~neg){
            int l=1,r=n-1;
            vec v=p[nxt(pos)]-p[pos];
            while(l<=r){
                int len=l+r>>1,mid=pos+len;
                if(mid>=n)mid-=n;
                if(!check(mid)){neg=mid;break;}
                int nx=nxt(mid);
                if(!check(nx)){neg=nx;break;}
                if(v.cross(p[nx]-p[mid])<eps)r=len-1;
                else l=len+1;
            }
            if(!~neg)return {0,n-1};
        }
        int L=pos,l=1,r=pos-neg-1;
        if(r<0)r+=n;
        while(l<=r){
            int len=l+r>>1,mid=pos-len;
            if(mid<0)mid+=n;
            if(check(mid))L=mid,l=len+1;
            else r=len-1;
        }
        int R=pos;l=1,r=neg-pos-1;
        if(r<0)r+=n;
        while(l<=r){
            int len=l+r>>1,mid=pos+len;
            if(mid>=n)mid-=n;
            if(check(mid))R=mid,l=len+1;
            else r=len-1;
        }
        R=nxt(R);
        if(L==R)return {0,n-1};
        return {L,R};
    }
    inline convex operator+(convex &rhs){
        if(!n||!rhs.n)return {};
        rotate(p.begin(),min_element(ALL(p),[&](const vec &u,const vec &v){return u.x!=v.x?u.x<v.x:u.y<v.y;}),p.end());
        rotate(rhs.p.begin(),min_element(ALL(rhs.p),[&](const vec &u,const vec &v){return u.x!=v.x?u.x<v.x:u.y<v.y;}),rhs.p.end());
        convex ret;
        ret.add(p[0]+rhs.p[0]);
        int i=0,j=0;
        vec nw=ret.p[0];
        while(i<n&&j<rhs.n){
            vec u=p[nxt(i)]-p[i],v=rhs.p[rhs.nxt(j)]-rhs.p[j];
            auto w=u.cross(v);
            if(w>0)++i,ret.add(nw=nw+u);
            else if(w<0)++j,ret.add(nw=nw+v);
            else ++i,++j,ret.add(nw=nw+u+v);
        }
        while(i<n)ret.add(nw=nw+p[nxt(i)]-p[i]),++i;
        while(j<rhs.n)ret.add(nw=nw+rhs.p[rhs.nxt(j)]-rhs.p[j]),++j;
        --ret.n,ret.p.qb();
        return ret;
    }
};

struct vec3D{
    double x,y,z;
    inline vec3D(double x=0.,double y=0.,double z=0.):x(x),y(y),z(z){}
    inline vec3D operator+(const vec3D &rhs){return {x+rhs.x,y+rhs.y,z+rhs.z};}
    inline vec3D operator-(const vec3D &rhs){return {x-rhs.x,y-rhs.y,z-rhs.z};}
    inline vec3D cross(const vec3D &rhs){return {y*rhs.z-z*rhs.y,z*rhs.x-x*rhs.z,x*rhs.y-y*rhs.x};}
    inline double dot(const vec3D &rhs){return x*rhs.x+y*rhs.y+z*rhs.z;}
    inline double operator*(const vec3D &rhs){return dot(rhs);}
    inline double norm(){return sqrt(x*x+y*y+z*z);}
    inline double projcoef(const vec3D &rhs){return dot(rhs)/dot(*this);}
    inline double projlen(const vec3D &rhs){return dot(rhs)/norm();}
};
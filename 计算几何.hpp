const double eps=1e-8,PI=acos(-1.);

struct vec2D{
    double x,y;
    inline vec2D(double x_=0,double y_=0):x(x_),y(y_){}
    inline vec2D operator+(const vec2D &rhs){return {x+rhs.x,y+rhs.y};}
    inline vec2D operator-(const vec2D &rhs){return {x-rhs.x,y-rhs.y};}
    inline vec2D operator*(double k){return {x*k,y*k};}
    inline double cross(const vec2D &rhs)const{return x*rhs.y-y*rhs.x;}
    inline double dot(const vec2D &rhs){return x*rhs.x+y*rhs.y;}
    inline double operator*(const vec2D &rhs){return dot(rhs);}
    inline bool operator==(const vec2D &rhs){return max(fabs(x-rhs.x),fabs(y-rhs.y))<eps;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2D &rhs){return dot(rhs)*dot(rhs)==dot(*this)*(rhs*rhs);}
    inline bool coln(const vec2D &rhs)const{return fabs(cross(rhs))<eps;}
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
    inline vec2d(ll x_=0,ll y_=0):x(x_),y(y_){}
    inline vec2d operator+(const vec2d &rhs){return {x+rhs.x,y+rhs.y};}
    inline vec2d operator-(const vec2d &rhs){return {x-rhs.x,y-rhs.y};}
    inline vec2d operator*(ll k){return {x*k,y*k};}
    inline ll cross(const vec2d &rhs)const{return x*rhs.y-y*rhs.x;}
    inline ll dot(const vec2d &rhs){return x*rhs.x+y*rhs.y;}
    inline ll operator*(const vec2d &rhs){return dot(rhs);}
    inline bool operator==(const vec2d &rhs){return x==rhs.x&&y==rhs.y;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2D &rhs){return dot(rhs)*dot(rhs)==dot(*this)*(rhs*rhs);}
    inline bool coln(const vec2d &rhs){return !cross(rhs);}
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
inline vec2D intersect(const vec &p1,const vec &v1,const vec &p2,const vec &v2){ // (p1, v1) and (p2, v2) are lines
    auto denom=v1.cross(v2);
    assert(abs(denom)>eps);
    double t=1.*(p2-p1).cross(v2)/denom;
    return vec2D(p1.x+v1.x*t,p1.y+v1.y*t);
}

template<class vec>
inline void angle_sort(V<vec> &v){
    sort(ALL(v),[&](const vec &x,const vec &y){
        return x.quad()!=y.quad()?
               x.quad()<y.quad():
               x.dir(y)==1; // eps may break sort, use atan2 instead
    });
}

template<class vec>
struct convex{
    V<vec>p;
    inline int nxt(int k){return k+1==p.size()?0:k+1;}
    inline int pre(int k){return k?k-1:p.size()-1;} // p should not be empty
    inline void insert(const vec &k){p.pb(k);}
    inline void init(bool coln=false){
        sort(ALL(p),[&](const vec &u,const vec &v){return u.x!=v.x?u.x<v.x:u.y<v.y;});
        p.erase(unique(ALL(p)),p.end());
        if(p.size()<3)return;
        V<vec>hi,lo;
        auto bad=[&](const V<vec> &v,const vec &k){
            if(v.size()<2)return false;
            auto cross=[&](const vec &a,const vec &b,const vec &c){
                // (b-a).cross(c-a)
                return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
            };
            auto c=cross(v[v.size()-2],v.back(),k);
            return coln?c<0:c<=0;
        };
        For(i,p.size()){
            while(bad(lo,p[i]))lo.qb();
            lo.pb(p[i]);
        }
        Rep(i,p.size()){
            while(bad(hi,p[i]))hi.qb();
            hi.pb(p[i]);
        }
        lo.qb(),hi.qb();
        p.swap(lo),p.insert(p.end(),ALL(hi));
    }
    inline convex(const V<vec> &p):p(p){init();}
    inline double area(){
        if(p.size()<3)return 0.;
        double ret=0.;
        For(i,p.size()){
            int j=nxt(i);
            if(j==p.size())j=0;
            ret+=p[i].cross(p[j]);
        }
        return ret*.5;
    }
    inline double peri(){
        if(p.size()<2)return 0.;
        double ret=0.;
        For(i,p.size()){
            int j=nxt(i);
            if(j==p.size())j=0;
            ret+=p[i].dist(p[j]);
        }
        return ret;
    }
    inline double diam(){
        if(p.size()<2)return 0.;
        if(p.size()==2)return p[0].dist(p[1]);
        auto cross=[&](const vec &a,const vec &b,const vec &c){
            // (b-a).cross(c-a)
            return (b.x-a.x)*(c.y-a.y)-(b.y-a.y)*(c.x-a.x);
        };
        int j=1;
        double ret=0;
        For(i,p.size()){
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
};

struct vec3D{
    double x,y,z;
    inline vec3D(double x_=0.,double y_=0.,double z_=0.):x(x_),y(y_),z(z_){}
    inline vec3D operator+(const vec3D &rhs){return {x+rhs.x,y+rhs.y,z+rhs.z};}
    inline vec3D operator-(const vec3D &rhs){return {x-rhs.x,y-rhs.y,z-rhs.z};}
    inline vec3D cross(const vec3D &rhs){return {y*rhs.z-z*rhs.y,z*rhs.x-x*rhs.z,x*rhs.y-y*rhs.x};}
    inline double dot(const vec3D &rhs){return x*rhs.x+y*rhs.y+z*rhs.z;}
    inline double operator*(const vec3D &rhs){return dot(rhs);}
    inline double norm(){return sqrt(x*x+y*y+z*z);}
    inline double projcoef(const vec3D &rhs){return dot(rhs)/dot(*this);}
    inline double projlen(const vec3D &rhs){return dot(rhs)/norm();}
};
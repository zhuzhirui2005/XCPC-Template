const double eps=1e-8;
struct vec2D{
    double x,y;
    inline vec2D(double x_=0,double y_=0):x(x_),y(y_){}
    inline vec2D operator+(const vec2D &rhs)const{return {x+rhs.x,y+rhs.y};}
    inline vec2D operator-(const vec2D &rhs)const{return {x-rhs.x,y-rhs.y};}
    inline double cross(const vec2D &rhs){return x*rhs.y-y*rhs.x;}
    inline double operator*(const vec2D &rhs)const{return x*rhs.x+y*rhs.y;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2D &rhs){return (*this*rhs)*(*this*rhs)==(*this**this)*(rhs*rhs);}
    inline bool coln(const vec2D &rhs){return fabs(cross(rhs))<eps;}
    inline int dir(const vec2D &rhs){return !coln(rhs)?cross(rhs)>=eps?1:-1:0;}
    inline double norm(){return sqrt(x*x+y*y);}
    inline double proj(const vec2D &rhs){return 1.*(*this*rhs)/(*this**this);}
    inline vec2D rot(double theta){
        double c=cos(theta),s=sin(theta);
        return {x*c+y*s,y*c-x*s};
    }
    inline int quad(){
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
    inline vec2d operator+(const vec2d &rhs)const{return {x+rhs.x,y+rhs.y};}
    inline vec2d operator-(const vec2d &rhs)const{return {x-rhs.x,y-rhs.y};}
    inline ll cross(const vec2d &rhs){return x*rhs.y-y*rhs.x;}
    inline ll operator*(const vec2d &rhs)const{return x*rhs.x+y*rhs.y;}
    // unsafe since overflow, use !cross() instead
    // inline bool coln(const vec2d &rhs){return (*this*rhs)*(*this*rhs)==(*this**this)*(rhs*rhs);}
    inline bool coln(const vec2d &rhs){return !cross(rhs);}
    inline int dir(const vec2d &rhs){return cross(rhs)?cross(rhs)>0?1:-1:0;}
    inline double norm(){return sqrt(x*x+y*y);}
    inline double proj(const vec2d &rhs){return 1.*(*this*rhs)/(*this**this);}
    inline vec2D rot(double theta){
        double c=cos(theta),s=sin(theta);
        return {x*c+y*s,y*c-x*s};
    }
    inline int quad(){
        if(x>0&&y>=0)return 1;
        if(x<=0&&y>0)return 2;
        if(x<0&&y<=0)return 3;
        if(x>=0&&y<0)return 4;
        return 0;
    }
};
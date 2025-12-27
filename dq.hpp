template<class T>
struct dq{
    T *a;
    int cap,h,n;
    inline dq():a(nullptr),cap(0),h(0),n(0){}
    inline ~dq(){delete[] a;}
    inline dq(dq &&k)noexcept:a(k.a),cap(k.cap),h(k.h),n(k.n){k.a=nullptr,k.cap=k.h=k.n=0;}
    inline dq &operator=(dq &&k)noexcept{swap(k);return *this;}
    inline void extend(){
        T *b=new T[cap?cap<<1:1];
        For(i,n)b[i]=a[(h+i)&(cap-1)];
        delete[] a;
        a=b,cap=cap?cap<<1:1,h=0;
    }
    inline void push_back(const T &k){if(n==cap)extend();a[(h+n++)&(cap-1)]=k;}
    inline void push_front(const T &k){if(n==cap)extend();a[h=(h-1)&(cap-1)]=k,++n;}
    inline void pop_back(){--n;}
    inline void pop_front(){h=(h+1)&(cap-1),--n;}
    inline const T &operator[](int k)const{return a[(h+k)&(cap-1)];}
    inline T &operator[](int k){return a[(h+k)&(cap-1)];}
    inline int size()const{return n;}
    inline void swap(dq &k){::swap(a,k.a),::swap(cap,k.cap),::swap(h,k.h),::swap(n,k.n);}
};
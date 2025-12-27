template<class T,class U=less<T>>
struct delpq{
    const U cmp=U();
    priority_queue<T,V<T>,U>q1,q2;
    inline delpq():q1(cmp),q2(cmp){}
    inline void push(const T &x){q1.push(x);}
    inline void pop(const T &x){q2.push(x);}
    inline void pop(){
        while(q2.size()&&q1.top()==q2.top())q1.pop(),q2.pop();
        assert(q1.size());
        q1.pop();
    }
    inline T top(){
        while(q2.size()&&q1.top()==q2.top())q1.pop(),q2.pop();
        assert(q1.size());
        return q1.top();
    }
    inline bool empty(){return q1.size()==q2.size();}
    inline int size(){assert(q1.size()>=q2.size());return q1.size()-q2.size();}
};

template<class T>
struct kpq{
    int k;
    delpq<T,greater<T>>s1;
    delpq<T>s2;
    ll sum; // 暂时默认对堆中前 k 大的数求和
    inline kpq(int k=0):k(k),sum(0){static_assert((is_same_v<T,int>)||(is_same_v<T,ll>));}
    inline void push(const T &x){
        if(s1.size()<k)s1.push(x),sum+=x;
        else{
            if(s1.size()&&s1.top()<x)s2.push(s1.top()),sum-=s1.top(),s1.pop(),s1.push(x),sum+=x;
            else s2.push(x);
        }
    }
    inline void pop(const T &x){
        if(s1.empty()||x<s1.top())s2.pop(x);
        else{
            sum-=x,s1.pop(x);
            if(s1.size()<k&&s2.size())s1.push(s2.top()),sum+=s2.top(),s2.pop();
        }
    }
    inline void reset(int kk){
        k=kk;
        while(s1.size()<k&&s2.size())s1.push(s2.top()),sum+=s2.top(),s2.pop();
        while(s1.size()>k)s2.push(s1.top()),sum-=s1.top(),s1.pop();
    }
};
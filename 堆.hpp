template<class T,class U=less<T>>
struct delpq{
    priority_queue<T,V<T>,U>q1,q2;
    inline delpq(){}
    inline delpq(const U &func){priority_queue<T,V<T>,U>(func).swap(q1),priority_queue<T,V<T>,U>(func).swap(q2);}
    inline void push(const T &x){q1.push(x);}
    inline void pop(const T &x){q2.push(x);}
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
    multiset<T>s1,s2;
    ll sum;
    inline kpq(int _k=0):k(_k),sum(0){}
    inline void insert(const T &x){
        if(s1.size()<k)s1.insert(x),sum+=x;
        else{
            if(x>*s1.begin())s2.insert(*s1.begin()),sum-=*s1.begin(),s1.erase(s1.begin()),s1.insert(x),sum+=x;
            else s2.insert(x);
        }
    }
    inline void erase(const T &x){
        if(s1.size()&&x<*s1.begin()){
            auto it=s2.find(x);
            assert(it!=s2.end());
            s2.erase(it);
        }
        else{
            auto it=s1.find(x);
            assert(it!=s1.end());
            s1.erase(it);
            sum-=x;
            if(s1.size()<k&&s2.size())s1.insert(*s2.rbegin()),sum+=*s2.rbegin(),s2.erase(prev(s2.end()));
        }
    }
};
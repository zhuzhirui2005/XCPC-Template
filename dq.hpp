template<class T>
struct dq{
    int hd;
    V<T>q;
    inline dq(){hd=0;}
    inline T front(int k=0){assert(hd+k<q.size());return q[hd+k];}
    inline T back(int k=0){assert(hd+k<q.size());return q[q.size()-1-k];}
    inline int size(){return q.size()-hd;}
    inline void clear(){hd=0,V<T>().swap(q);}
    inline void push(const T &v){q.pb(v);}
    inline void pop_back(){q.qb();}
    inline void pop_front(){assert(hd<q.size());++hd;}
};
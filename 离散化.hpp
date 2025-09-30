template<class T>
struct disc{
	// 0-indexed
	vector<T>d;
	inline disc(){}
	inline void insert(const T &x){d.pb(x);}
	inline void insert(const V<T> &v){d.insert(d.end(),ALL(v));}
	inline void init(){sort(ALL(d));d.erase(unique(ALL(d)),d.end());}
	inline disc(const vector<T> &v){d=move(v);init();}
	inline int query(const T &x){return lower_bound(ALL(d),x)-d.begin();}
	inline int size(){return d.size();}
};

template<class T,class container>
struct genID{
    int cnt;
    container id;
    inline genID():cnt(0){}
    inline bool count(const T &ele){return id.count(ele);}
    inline int get_id(const T &ele){
        auto it=id.emplace(ele,-1);
        if(it.se)it.fi->se=cnt++;
        return it.fi->se;
    }
};
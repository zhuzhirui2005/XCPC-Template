#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

template<class T>
struct rbt{
	typedef pair<T,int> pti;
	int cnt;
	typedef tree<pti,null_type,less<pti>,rb_tree_tag,tree_order_statistics_node_update> rbt_t;
	rbt_t t;
	inline rbt(){cnt=0;}
	inline void clear(){cnt=0,rbt_t().swap(t);}
	inline typename rbt_t::iterator begin(){return t.begin();}
	inline typename rbt_t::iterator end(){return t.end();}
	inline void insert(const T &x){t.insert({x,cnt++});}
	inline typename rbt_t::iterator find(const T &x){return t.lower_bound({x,0});}
	inline void erase(const T &x){t.erase(find(x));}
	inline T pre(const T &x){
		auto it=find(x);
		assert(it!=begin());
		return prev(it)->fi;
	}
	inline T nxt(const T &x){
		auto it=find(x+1);
		assert(it!=end());
		return it->fi;
	}
	// all 0-indexed
	inline int rk(const T &x){return t.order_of_key({x,0});}
	inline T at(unsigned x){return t.find_by_order(x)->fi;}
};

#include <ext/pb_ds/priority_queue.hpp>
inline V<ll> dijkstra(int n,int s,const V<V<pii>> &to){
	assert(0<=n),assert(0<=s),assert(s<n),assert(to.size()<=n);
	for(const V<pii> &i:to)
        for(const pii &j:i)
            assert(0<=min(j.fi,j.se)),assert(j.fi<n);
    V<ll>dis(n,infl);
    dis[s]=0;
    __gnu_pbds::priority_queue<pli,greater<pli>,pairing_heap_tag>q;
    V<decltype(q)::point_iterator>it(n);
    it[s]=q.push({0,s});
    while(q.size()){
        int p=q.top().se;q.pop();
        for(const pii &i:to[p])
            if(ckmin(dis[i.fi],dis[p]+i.se)){
				if(it[i.fi]!=NULL)q.modify(it[i.fi],{dis[i.fi],i.fi});
                else it[i.fi]=q.push({dis[i.fi],i.fi});
            }
    }
    for(ll &i:dis)if(i==infl)i=-1;
    return dis;
}
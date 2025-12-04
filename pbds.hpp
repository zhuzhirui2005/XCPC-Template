#include <ext/pb_ds/tree_policy.hpp>
#include <ext/pb_ds/assoc_container.hpp>
using namespace __gnu_pbds;

template<class NCIT,class NIT,class CMP,class ALC>
struct custom_node_update:detail::branch_policy<NCIT,NIT,ALC>{
    typedef detail::branch_policy<NCIT,NIT,ALC> base_type;
    using CIT=typename NCIT::value_type;
    using size_type=typename ALC::size_type;
    typedef typename base_type::key_type key_type;
    using metadata_type=pli;
    virtual NCIT node_begin()const=0;
    virtual NCIT node_end()const=0;
    virtual CMP& get_cmp_fn()=0;
    inline void operator()(NIT it,NCIT ed){
        NIT l=it.get_l_child(),r=it.get_r_child();
        metadata_type nw={(*it)->fi,1};
        if(l!=ed){
            metadata_type L=l.get_metadata();
            nw.fi+=L.fi,nw.se+=L.se;
        }
        if(r!=ed){
            metadata_type R=r.get_metadata();
            nw.fi+=R.fi,nw.se+=R.se;
        }
        const_cast<metadata_type&>(it.get_metadata())=nw;
    }
    inline ll prefix_sum_by_order(size_type k){
        if(k<0)return 0;
        ++k;
        NCIT it=node_begin(),ed=node_end();
        ll ret=0;
        while(it!=ed&&k){
            NCIT l=it.get_l_child();
            if(l!=ed){
                metadata_type L=l.get_metadata();
                if(k<=L.se){
                    it=l;
                    continue;
                }
                ret+=L.fi,k-=L.se;
            }
            ret+=(*it)->fi;
            if(!--k)break;
            it=it.get_r_child();
        }
        return ret;
    }
    CIT find_by_order(size_type k)const{
        ++k;
        NCIT it=node_begin(),ed=node_end();
        while(it!=ed){
            NCIT l=it.get_l_child();
            if(l!=ed){
                metadata_type L=l.get_metadata();
                if(k<=L.se){
                    it=l;
                    continue;
                }
                k-=L.se;
            }
            if(!--k)return *it;
            it=it.get_r_child();
        }
        return base_type::end_iterator();
    }
    size_type order_of_key(const key_type &x){
        NCIT it=node_begin(),ed=node_end();
        CMP &cmp=get_cmp_fn();
        size_type ret=0;
        while(it!=ed){
            NCIT l=it.get_l_child();
            if(cmp(x,**it)){
                it=l;
                continue;
            }
            if(l!=ed){
                metadata_type L=l.get_metadata();
                ret+=L.se;
            }
            if(cmp(**it,x))it=it.get_r_child(),++ret;
            else break;
        }
        return ret;
    }
};

template<class T,template<class,class,class,class>class node_update=tree_order_statistics_node_update>
struct rbt{
	typedef pair<T,int> pti;
	int cnt;
	typedef tree<pti,null_type,less<pti>,rb_tree_tag,node_update> rbt_t;
	rbt_t t;
	inline rbt():cnt(0){}
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
				if(it[i.fi]!=nullptr)q.modify(it[i.fi],{dis[i.fi],i.fi});
                else it[i.fi]=q.push({dis[i.fi],i.fi});
            }
    }
    for(ll &i:dis)if(i==infl)i=-1;
    return dis;
}
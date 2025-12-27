namespace _Uncategorized_1{
    // 三维偏序
    template<class T>
    inline V<pii> solve(int n,int m,const V<T> &x,const V<int> &y){
        // 对每个 i 求左右两边第一个 x_j >= x_i 且 y_j >= y_i 的 j，y 需要被离散化
        assert(x.size()==n),assert(y.size()==n);
        for(int i:y)assert(0<=i),assert(i<m);
        int cnt=0;
        V<int>c(m+1),id(n),ver(m+1);
        iota(ALL(id),0);
        V<pii>ret(n,{-1,n});
        auto cdq=[&](auto &&self,int l,int r)->void{
            if(l==r)return;
            int mid=l+r>>1;
            self(self,l,mid),self(self,mid+1,r);
            ++cnt;
            for(int i=l,j=mid+1;j<=r;++j){
                while(i<=mid&&x[id[i]]>=x[id[j]]){
                    for(int k=y[id[i]]+1;k;k&=k-1){
                        if(ver[k]<cnt)c[k]=id[i],ver[k]=cnt;
                        else ckmax(c[k],id[i]);
                    }
                    ++i;
                }
                for(int k=y[id[j]]+1;k<=m;k+=k&-k)if(ver[k]==cnt)ckmax(ret[id[j]].fi,c[k]);
            }
            ++cnt;
            for(int i=l,j=mid+1;i<=mid;++i){
                while(j<=r&&x[id[j]]>=x[id[i]]){
                    for(int k=y[id[j]]+1;k;k&=k-1){
                        if(ver[k]<cnt)c[k]=id[j],ver[k]=cnt;
                        else ckmin(c[k],id[j]);
                    }
                    ++j;
                }
                for(int k=y[id[i]]+1;k<=m;k+=k&-k)if(ver[k]==cnt)ckmin(ret[id[i]].se,c[k]);
            }
            inplace_merge(id.begin()+l,id.begin()+mid+1,id.begin()+r+1,[&](int u,int v){return x[u]>x[v];});
        };
        cdq(cdq,0,n-1);
        return ret;
    }
}

namespace _Uncategorized_2{
    // 贪心
    template<class T,class U>
    inline bool solve(int n,T m,const V<U> &a,const V<U> &b){
        // n 个任务，初始有 m 的资金，第 i 个任务需要资金 >= a_i 才能完成，完成后获得 b_i 的资金，判断能否全部完成
        static_assert((is_unsigned_v<T>)==(is_unsigned_v<U>));
        assert(a.size()==n),assert(b.size()==n);
        V<int>id(n);
        iota(ALL(id),0);
        {
            auto it=partition(ALL(id),[&](int k){return b[k]>=0;});
            sort(id.begin(),it,[&](int x,int y){return a[x]<a[y];});
            sort(it,id.end(),[&](int x,int y){return a[x]+b[x]>a[y]+b[y];});
        }
        for(int i:id){
            if(m<a[i])return false;
            m+=b[i];
        }
        return true;
    }
}

namespace _Uncategorized_3{
}

namespace _Uncategorized_4{
    // 扫描线
    inline V<array<int,3>> solve(const V<int> &a){
        // 求 a 的所有极小 mex 区间，mex 为未出现过的最小非负整数
        int n=a.size();
        V<int>lst(n,-1);
        V<array<int,3>>ret;
        V<int>t(n<<2,-1);
        auto modify=[&](auto &&self,int p,int l,int r,int k,int v)->int{
            int ret;
            if(l==r)ret=t[p],t[p]=v;
            else{
                int mid=l+r>>1;
                ret=k<=mid?self(self,p<<1,l,mid,k,v):self(self,p<<1|1,mid+1,r,k,v);
                t[p]=min(t[p<<1],t[p<<1|1]);
            }
            return ret;
        };
        auto query=[&](auto &&self,int p,int l,int r,int k)->int{
            if(r<=k)return t[p];
            int mid=l+r>>1;
            return k<=mid?self(self,p<<1,l,mid,k):min(t[p<<1],self(self,p<<1|1,mid+1,r,k));
        };
        auto find=[&](auto &&self,int p,int l,int r,int k)->int{
            if(l==r)return l;
            int mid=l+r>>1;
            return t[p<<1]<k?self(self,p<<1,l,mid,k):self(self,p<<1|1,mid+1,r,k);
        };
        For(i,n){
            assert(a[i]>=0);
            if(a[i])ret.pb({i,i,0});
            if(a[i]>=n)continue;
            int j=modify(modify,1,0,n-1,a[i],i),k=query(query,1,0,n-1,a[i]);
            lst[a[i]]=i;
            while(j<k){
                int l=t[1]<k?find(find,1,0,n-1,k):n;
                ret.pb({k,i,l});
                if(l==n)break;
                k=lst[l];
            }
        }
        return ret;
    }
}

namespace _Uncategorized_5{
    // 长剖 DP
    template<class T,int k=2>
    inline T solve(int K,const V<V<int>> &to){
        // 在 n 个点的树中选 k 个点使得 k 个点两两距离均相同的方案数
        if constexpr(k<11){
            if(k<K)return solve<T,k+1>(K,to);
        }
        assert(k==K);
        int n=to.size();
        if(k==2)return (1ll*n*(n-1)>>1)%mod;
        V<int>mx(n),son(n,-1);
        auto dfs1=[&](auto &&self,int p,int fa)->void{
            for(int i:to[p])if(i!=fa){
                self(self,i,p);
                if(ckmax(mx[p],mx[i]))son[p]=i;
            }
            ++mx[p];
        };
        dfs1(dfs1,0,-1);
        T ret=0;
        auto dfs2=[&](auto &&self,int p,int fa)->array<dq<T>,k-1>{
            array<dq<T>,k-1>f;
            if(~son[p]){
                array<dq<T>,k-1>g=self(self,son[p],p);
                if(g[k-2].size()){
                    g[k-2].qf();
                    if(g[k-2].size())ret+=g[k-2][0];
                }
                f[0].swap(g[0]),f[k-2].swap(g[k-2]);
            }
            f[0].pf(1);
            assert(f[0].size()==mx[p]);
            for(int i:to[p])if(i!=fa&&i!=son[p]){
                array<dq<T>,k-1>g=self(self,i,p);
                g[0].pf(0);
                if(g[k-2].size())g[k-2].qf();
                For(j,min(f[0].size(),g[k-2].size()))ret+=f[0][j]*g[k-2][j];
                For(j,min(f[k-2].size(),g[0].size()))ret+=f[k-2][j]*g[0][j];
                if(f[k-2].size()<g[k-2].size())f[k-2].swap(g[k-2]);
                For(j,g[k-2].size())f[k-2][j]+=g[k-2][j];
                Rep(j,k-2)For(l,min(f[j].size(),g[0].size())){
                    if(l<f[j+1].size())f[j+1][l]+=f[j][l]*g[0][l];
                    else f[j+1].pb(f[j][l]*g[0][l]);
                }
                if(f[0].size()<g[0].size())f[0].swap(g[0]);
                For(j,g[0].size())f[0][j]+=g[0][j];
            }
            return f;
        };
        dfs2(dfs2,0,-1);
        return ret;
    }
}
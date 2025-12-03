inline V<pii> order3d(int n,int m,const V<int> &x,const V<int> &y){
    // 对每个 i 求左右两边第一个 x_j >= x_i 且 y_j >= y_i 的 j，y 需要被离散化
    assert(x.size()==n),assert(y.size()==n);
    for(int i:y)assert(0<=i),assert(i<m);
    int cnt=0;
    V<int>c(m+1),id(n),ver(m+1);
    iota(ALL(id),0);
    V<pii>ret(n,{-1,n});
    auto solve=[&](auto &&self,int l,int r)->void{
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
    solve(solve,0,n-1);
    return ret;
}
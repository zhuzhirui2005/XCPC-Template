struct dsu{
    int n;
    V<int>fa;
    inline dsu(int n=0):n(n),fa(n,-1){}
    int find(int k){return fa[k]<0?k:fa[k]=find(fa[k]);}
    inline bool merge(int x,int y){
        x=find(x),y=find(y);
        if(x!=y)fa[x]+=fa[y],fa[y]=x;
        return x!=y;
    }
    inline bool same(int x,int y){return find(x)==find(y);}
    inline int size(int k){return -fa[find(k)];}
};

struct range_dsu{
    V<V<int>>fa;
    int lg,n;
    inline range_dsu(int n=0):n(n){
        if(n){
            fa.resize(lg=__lg(n)+1);
            For(i,lg)fa[i].resize(n-(1<<i)+1,-1);
        }
        else lg=0;
    }
    int find(int d,int k){return fa[d][k]<0?k:fa[d][k]=find(d,fa[d][k]);}
    inline void merge(int d,int x,int y){
        x=find(d,x),y=find(d,y);
        if(x>y)swap(x,y);
        if(x!=y)fa[d][x]+=fa[d][y],fa[d][y]=x;
    }
    inline void merge(int x1,int x2,int y1,int y2){
        assert(x2-x1==y2-y1);
        Rep(i,lg)if(x1+(1<<i)-1<=x2){
            merge(i,x1,y1);
            x1+=1<<i,y1+=1<<i;
        }
    }
    inline void init(){
        REP(i,1,lg)For(j,n-(1<<i)+1){
            int k=find(i,j);
            merge(i-1,j,k),merge(i-1,j+(1<<i-1),k+(1<<i-1));
        }
    }
    int find(int k){return fa[0][k]<0?k:fa[0][k]=find(fa[0][k]);}
    inline bool same(int x,int y){return find(x)==find(y);}
    inline int size(int k){return -fa[0][find(k)];}
};

struct pal_dsu{
    range_dsu d;
    int n;
    inline pal_dsu(int n=0):d(n<<1),n(n){For(i,n)d.merge(i,i,(n<<1)-1-i,(n<<1)-1-i);}
    inline void merge(int x1,int x2,int y1,int y2){d.merge(x1,x2,(n<<1)-1-y2,(n<<1)-1-y1);}
    inline void init(){d.init();}
    int find(int k){return d.find(k);}
    inline bool same(int x,int y){return d.same(x,y);}
    inline int size(int k){return d.size(k)>>1;}
};
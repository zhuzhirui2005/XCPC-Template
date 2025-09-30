template<int n>
struct KDT{
    array<int,n>mn,mx,t;
    KDT *ls,*rs;
    inline KDT(const array<int,n> &arr={}):mn(arr),mx(arr),t(arr),ls(NULL),rs(NULL){}
    inline void build(V<array<int,n>> &a,int l,int r,int d){
        int mid=l+r>>1;
        nth_element(a.begin()+l,a.begin()+mid,a.begin()+r+1,[&](const auto &x,const auto &y){return x[d]<y[d];});
        *this=KDT(a[mid]);
        if(l<mid){
            ls=new KDT();
            ls->build(a,l,mid-1,d+1<n?d+1:0);
            For(i,n)ckmin(mn[i],ls->mn[i]),ckmax(mx[i],ls->mx[i]);
        }
        if(r>mid){
            rs=new KDT();
            rs->build(a,mid+1,r,d+1<n?d+1:0);
            For(i,n)ckmin(mn[i],rs->mn[i]),ckmax(mx[i],rs->mx[i]);
        }
    }
};
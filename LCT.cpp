struct LCT{
    int n;
    V<int>fa,siz;
    V<bool>rev;
    V<mi>sum,val;
    V<array<int,2>>son;
    V<pair<mi,mi>>tag;
    inline LCT(int n):fa(n,-1),n(n),rev(n),siz(n,1),son(n,{-1,-1}),sum(n,1),tag(n,{1,0}),val(n,1){}
    inline bool is_rt(int p){return !~fa[p]||(son[fa[p]][0]!=p&&son[fa[p]][1]!=p);}
    inline void push_up(int p){
        siz[p]=1,sum[p]=val[p];
        if(~son[p][0])siz[p]+=siz[son[p][0]],sum[p]+=sum[son[p][0]];
        if(~son[p][1])siz[p]+=siz[son[p][1]],sum[p]+=sum[son[p][1]];
    }
    inline void apply_rev(int p){rev[p]=!rev[p],swap(son[p][0],son[p][1]);}
    inline void apply_tag(int p,const pair<mi,mi> &v){
        sum[p]=v.fi*sum[p]+siz[p]*v.se;
        tag[p]={v.fi*tag[p].fi,v.fi*tag[p].se+v.se};
        val[p]=v.fi*val[p]+v.se;
    }
    inline void push_down(int p){
        if(rev[p]){
            if(~son[p][0])apply_rev(son[p][0]);
            if(~son[p][1])apply_rev(son[p][1]);
            rev[p]=false;
        }
        if(tag[p].fi!=1||tag[p].se.val){
            if(~son[p][0])apply_tag(son[p][0],tag[p]);
            if(~son[p][1])apply_tag(son[p][1],tag[p]);
            tag[p]={1,0};
        }
    }
    inline void rotate(int p){
        int f=fa[p],g=fa[f];
        bool r=son[f][1]==p;
        int &q=son[p][!r];
        if(!is_rt(f))son[g][son[g][1]==f]=p;
        son[f][r]=q;
        if(~q)fa[q]=f;
        fa[q=f]=p,fa[p]=g;
        push_up(f),push_up(p);
    }
    inline void splay(int p){
        auto dfs=[&](auto &&self,int p)->void{
            if(!is_rt(p))self(self,fa[p]);
            push_down(p);
        };
        dfs(dfs,p);
        while(!is_rt(p)){
            int f=fa[p],g=fa[f];
            if(!is_rt(f))rotate((son[f][0]==p)==(son[g][0]==f)?f:p);
            rotate(p);
        }
    }
    inline void access(int p){
        int s=-1;
        while(~p){
            splay(p);
            son[p][1]=s,push_up(p);
            s=p,p=fa[p];
        }
    }
    inline void modify(int p,mi v){
        access(p),splay(p);
        sum[p]+=v-val[p],val[p]=v;
        push_up(p);
    }
    inline void make_root(int p){access(p),splay(p),apply_rev(p);}
    inline void split(int x,int y){make_root(x),access(y),splay(y);}
    inline void add(int x,int y,const pair<mi,mi> &v){split(x,y),apply_tag(y,v);}
    inline mi qry(int x,int y){split(x,y);return sum[y];}
    inline int find_root(int p){
        access(p),splay(p);
        while(~son[p][0])push_down(p),p=son[p][0];
        splay(p);
        return p;
    }
    inline bool link(int x,int y){
        make_root(x);
        if(find_root(y)==x)return false;
        fa[x]=y;
        return true;
    }
    inline bool cut(int x,int y){
        split(x,y);
        if(~son[x][1]||son[y][0]!=x)return false;
        fa[x]=son[y][0]=-1;
        push_up(y);
        return true;
    }
};
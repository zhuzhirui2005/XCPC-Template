#include<bits/stdc++.h>
using namespace std;
int main(){
    long long k,p,q;
    scanf("%lld%lld%lld",&p,&q,&k);
    int ans=0;
    for(int i=1;i<=q;++i){
        int j=i/gcd<int>(p,i);
        for(int l=2;l<=k;++l)if(k%l==0)while(j%l==0)j/=l;
        if(j>1)++ans;
    }
    printf("%d\n",ans);
}
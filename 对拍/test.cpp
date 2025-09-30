#include<bits/stdc++.h>
using namespace std;
int main(){
    long long k,p,q;
    scanf("%lld%lld%lld",&p,&q,&k);
    vector<long long>a{1};
    long long t=k;
    for(int i=2,j=sqrtl(k);i<=j;++i)if(k%i==0){
        do k/=i;while(k%i==0);
        int tmp=a.size();
        for(int l=0;l<tmp&&a[l]<=q/i;++l){
            a.push_back(a[l]*i);
            while(a.back()<=q/i)a.push_back(a.back()*i);
        }
        sort(a.begin(),a.end());
    }
    if(k>1){
        int tmp=a.size();
        for(int l=0;l<tmp&&a[l]<=q/k;++l){
            a.push_back(a[l]*k);
            while(a.back()<=q/k)a.push_back(a.back()*k);
        }
        sort(a.begin(),a.end());
    }
    long long ans=q;
    int l=a.size()-1;
    for(int i=1,j=sqrtl(p);i<=j&&l>=0;++i)if(p%i==0){
        while(l>=0&&i>q/a[l])--l;
        if(gcd(i,t)==1)ans-=l+1;
    }
    for(int i=sqrtl(p);i>0&&l>=0;--i)if(p%i==0){
        long long j=p/i;
        if(i!=j&&gcd(j,t)==1){
            while(l>=0&&j>q/a[l])--l;
            ans-=l+1;
        }
    }
    printf("%lld\n",ans);
}
/*
p/i 
*/
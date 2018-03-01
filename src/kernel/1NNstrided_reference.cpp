#include<stdio.h>
#include<cmath>

typedef double dt;
typedef long dti;
// this could implement pairwise later

// It could also be optimized to do sets of 8

//template<typename dt>
void initcov(dt *a, dt *cx, dt *mu, dt *invn, dt *mp, dti *mpi, int mindiag, int len, int subLen){
   dt m1 = mu[0];
   dt s1 = invn[0];
   dt maxmp = -1;
   dti maxmpi = -1;
   int upperlim = len - subLen + 1;
   for(int diag = mindiag; diag < upperlim; diag++){
      dt cov = 0;
      dt m2 = mu[diag];
      for(int offset = 0; offset < subLen; offset++){
         dt c1 = a[offset] - m1;
         dt c2 = a[diag+offset] - m2;
         cov = fma(c1,c2,cov);
      }
      cx[diag] = cov;
   }
}



static inline void solvetile2(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len) __attribute__ ((always_inline));

static inline void solvetile2(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len){
   for(int diag = dm; diag < dm + len; diag++){
      dt cov = cx[diag];
      for(int offset = om; offset < om + len; offset++){
         int k = diag+offset;
         cov = fma(df[offset],dx[k],cov);
         cov = fma(df[k],dx[offset],cov);
         dt corr = cov*s[offset];
         corr *= s[k];
         if(mp[offset] < corr){
            mp[offset] = corr;
            mpi[offset] = k;
         }
         if(mp[k] < corr){
            mp[k] = corr; 
            mpi[k] = offset;
         }
      }
      cx[diag] = cov;
   }
}




void solvempref(dt* __restrict__ a, dt*  __restrict__ cx, dt* __restrict__ mu, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int len, int diagmin, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);
   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
   for(int i = diagmin; i < len-sublen+1; i++){
      int jmax = len - sublen - i + 1;
      dt cov = cx[i];
      for(int j = 0; j < jmax; j++){
         int k = i+j;
         cov = fma(df[j],dx[k],cov); 
         cov = fma(df[k],dx[j],cov);
         dt cor = cov*s[j];
         cor *= s[k];
         if(mp[j] < cor){
            mp[j] = cor;
            mpi[j] = k;
         }
         if(mp[k] < cor){
            mp[k] = cor;
            mpi[k] = j;
         }
      }
      cx[i] = cov; 
   }
}




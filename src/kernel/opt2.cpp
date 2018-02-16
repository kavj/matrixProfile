#include<stdio.h>
#include<cmath>
#include<algorithm>
#include "opt1_utils.h"

typedef double dt;
typedef long dti;
// this could implement pairwise later

// It could also be optimized to do sets of 8

//template<typename dt>
static void initcov(const dt* __restrict__ a, dt* __restrict__ cx, const dt* __restrict__ mu, dt* invn, int len, int sublen);
static void initcov(const dt *a, dt *cx, const dt *mu, dt *invn, int minlag, int len, int subLen){
   a = (dt*)__builtin_assume_aligned(a,64);
   cx = (dt*)__builtin_assume_aligned(cx,64);
   mu = (dt*)__builtin_assume_aligned(mu,64);
   invn = (dt*)__builtin_assume_aligned(invn,64);

   dt m1 = mu[0];
   dt s1 = invn[0];
   dt maxmp = -1;
   dti maxmpi = -1;
   int upperlim = len - subLen + 1;
   for(int diag = minlag; diag < upperlim; diag++){
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


static inline void solvetile(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len) __attribute__ ((always_inline));

static inline void solvetile(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);

   for(int diag = dm; diag < dm + len; diag++){
      dt ca, cb, cc, cd, ce, cf, cg, ch;
      load_x8(ca, cb, cc, cd, ce, cf, cg, ch, cx, diag);
      for(int offset = om; offset < om + len; offset++){

         int k = diag+offset;

         dt ra = df[offset];
         fma_x8(ca, cb, cc, cd, ce, cf, cg, ch, ra, dx, offset); 

         ra = dx[offset];
         fma_x8(ca, cb, cc, cd, ce, cf, cg, ch, ra, df, offset);

         dt da, db, dc, dd, de, df, dg, dh;
         prod_oop_x8(da,db,dc,dd,de,df,dg,dh,ca,cb,cc,cd,ce,cf,cg,ch,s,k);

         ra = s[offset];
         prod_x8(da,db,dc,dd,de,df,dg,dh,ra);
         min_x8(da,db,dc,dd,de,df,dg,dh,mp,mpi,diag+offset,offset);
      }
      store_x8(ca, cb, cc, cd, ce, cf, cg, ch, cx, diag);
   }
}


static inline void solvefringe(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len) __attribute__ ((always_inline));

static inline void solvefringe(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int dmx, int omx){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);

   for(int diag = dm; diag < dmx; diag++){
      dt cov = cx[diag];
      int offmax = omx - diag;
      for(int offset = om; offset < offmax; offset++){
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


void solvempref(dt* __restrict__ a, dt*  __restrict__ cx, dt* __restrict__ mu, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int len, int minlag, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);

   initcov(a,cx,mu,s,minlag,len,sublen);
   int diagmax = len - sublen + 1;
   const int tsz = 128;
   for(int diag = minlag; diag < diagmax; diag += tsz){
      int cmp = diag + 2*tsz - 3;
      for(int offset = 0; offset < diagmax-diag; offset+=tsz){
         if(offset + cmp < diagmax){
            solvetile(cx,df,dx,s,mp,mpi,diag,offset,tsz);
         }
         else{
            int dmx = std::min(diagmax,diag+tsz);
            int omx = std::min(diagmax,offset+tsz);
            solvefringe(cx,df,dx,s,mp,mpi,diag,offset,dmx,omx);
         }
      }
   }
}



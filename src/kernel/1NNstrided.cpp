#include<stdio.h>
#include<cmath>
#include<algorithm>
#include "reg.h"
typedef double dt;
typedef long dti;
using namespace std;
// this could implement pairwise later

// It could also be optimized to do sets of 8


void initcov(dt *a, dt *cx, dt *mu, dt *invn, dt *mp, dti *mpi, int minlag, int len, int subLen){
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


// Todo: fix ordering so this gives the correct results
static inline void solvetile(dt* __restrict__cx, const dt* __restrict__ df, const dt* __restrict__dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len) __attribute__ ((always_inline));

static inline void solvetile(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);
   const int unroll = 10; 

   for(int diag = dm; diag < dm+len; diag+=unroll){
      reg_block<double> cov;
      for(int i = 0; i < unroll; i++){
         cov(i) = cx[diag+i];
      }
      for(int offset = om; offset < om+len; offset++){
         int k = offset+diag;
         dt ra = dx[offset];
         for(int i = 0; i < unroll; i++){
           cov(i) += ra*df[k+i];
         } 
         ra = df[offset];
         for(int i = 0; i < unroll; i++){
           cov(i) += ra*dx[k+i];
         } 
         reg_block<double> corr;
         for(int i = 0; i < unroll; i++){
            corr(i) = cov(i) * s[k+i];
         } 
         ra = s[offset];
         for(int i = 0; i < unroll; i++){
            corr(i) = cov(i) * ra;
         }
         for(int i = 0; i < unroll; i++){
            if(mp[k+i] < corr(i)){
            //   mpi[k+i] = offset;
               mp[k+i] = corr(i);
            }
         }
    //     reg_block<long> index;
         for(int i = 0; i < unroll; i+=2){
            if(corr(i) < corr(i+1)){
              // index(i) = k+i+1;
               corr(i) = corr(i+1);
            }
           /* else{
              index(i) = k+i;
            }*/
         }
         for(int i = 0; i < unroll; i+=4){
            if(corr(i) < corr(i+2)){
  //             index(i) = index(i+1);
               corr(i) = corr(i+2);
            }
         }
         if(corr(0) < corr(4)){
//            index(0) = index(4);
            corr(0) = corr(4);
            
         }
         if(mp[offset] < corr(0)){
          // mpi[offset] = index(0);
            mp[offset] = corr(0);
         }
      }
      for(int i = 0; i < unroll; i++){
         cx[diag+i] = cov(i);
      }
   }
}

void solvempref(dt *a, dt *cx, dt *mu, dt *df, dt *dx, dt *s, dt *mp, dti *mpi, int len, int diagmin, int sublen){
   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
   const int tilelen = 2051;
   for(int i = diagmin; i < len-sublen-tilelen; i+= tilelen){
      int jmax = len - sublen - i - tilelen + 1;
      for(int j = 1; j < jmax; j+= tilelen){
         solvetile(cx,df,dx,s,mp,mpi,i,j,tilelen);
      }
   }
}


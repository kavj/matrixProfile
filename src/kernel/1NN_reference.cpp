#include<stdio.h>
#include<cmath>

typedef double dt;
typedef long dti;
// this could implement pairwise later

// It could also be optimized to do sets of 8


void maxpearson_reference(dt*  __restrict__ cx, const dt* __restrict__ mu, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ invn, dt* __restrict__ mp, dti* __restrict__ mpi, int len, int diagmin, int sublen){
   for(int i = diagmin; i < len; i++){
      dt cov = cx[i];
      // if I don't update at initialization, I would need to do so here
      for(int j = i; j < len; j++){
         int k = j - i;
         cov = cov + df[j]*dx[k] + df[k]*dx[j];
         dt cor = cov*invn[j]*invn[k];
         if(mp[j] < cor){
            mp[j] = cor;
            mpi[j] = k;
         }
         if(mp[k] < cor){
            mp[k] = cor;
            mpi[k] = j;
         }
      }
   }
}


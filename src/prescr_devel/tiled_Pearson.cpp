#include<cstdio>
#include<algorithm>
#define klen 64 


static inline  __attribute__((always_inline)) void pauto_pearson_refkern (
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int offsetr,
   const int offsetc)
{
    cov =  (double*)__builtin_assume_aligned(cov,32);
    mp =   (double*)__builtin_assume_aligned(mp,32);
    mpi =  (int*)__builtin_assume_aligned(mpi,32);
    df =   (const double*)__builtin_assume_aligned(df,32);
    dg =   (const double*)__builtin_assume_aligned(dg,32);
    invn = (const double*)__builtin_assume_aligned(invn,32);

    
    for(int i = 0; i < klen; i++){
      for(int j = 0; j < klen; j++){
         cov[j] += dg[i]*df[i+j+offsetc];
      }
      double corr[klen];
      for(int j = 0; j < klen; j++){
         cov[j] += df[i]*dg[i+j+offsetc];
         corr[j] = cov[j]*invn[i]*invn[i+j+offsetc];
      }
      for(int j = 0; j < klen; j++){
         if(mp[j] < corr[j]){
            mp[j] = corr[j];
            mpi[j] = i+j+offsetr+offsetc;
         }
      }
      for(int j = 0; j < klen; j++){
         if(mp[i+j+offsetc] < corr[j]){
            mp[i+j+offsetc] = corr[j];
            mpi[i+j+offsetc] = j+offsetr;
         }
      }
   }
}


// annoying discrepancy between int and long long. May need to template these
// The problem is that int is generally ideal unless long long happens to match the size of a particular mask
static inline void pauto_pearson_edgekern(
   double*       __restrict__ cov, 
   double*       __restrict__ mp,  
   int*    __restrict__ mpi, 
   const double* __restrict__ df,  
   const double* __restrict__ dg, 
   const double* __restrict__ invn, 
   int offsetr, 
   int offsetc, 
   int bound)
{
    cov =  (double*)__builtin_assume_aligned(cov,32);
    mp =   (double*)__builtin_assume_aligned(mp,32);
    mpi =  (int*)__builtin_assume_aligned(mpi,32);
    df =   (const double*)__builtin_assume_aligned(df,32);
    dg =   (const double*)__builtin_assume_aligned(dg,32);
    invn = (const double*)__builtin_assume_aligned(invn,32);

 
   int dlim = std::min(klen,bound);
   for(int i = 0; i < dlim; i++){
      double c = cov[i];
      // make tighter bound later
      for(int j = i; j < bound-i; j++){
         c += df[j]*dg[j+offsetc];
         c += df[j+offsetc]*dg[j];
         double corr = c*invn[i]*invn[i+j];
         if(corr > mp[i]){
            mp[j] = corr;
            mpi[j] = j+offsetr+offsetc;
         }
         if(corr > mp[i+j+offsetc]){
            mp[j+offsetc] = corr;
            mpi[j+offsetc] = j+offsetr;
         }
      }
      cov[i] = c;
   }
}

//Todo: make preambles compiler agnostic


// offsetr indicates the amount which we have already offset all pointers other than cov from their base values
// offsetc indicates the amount which we have offset cov from its base value
// Since we require 2 values from mp,mpi,df,dg,invn per update, we apply offsetc as needed to index calculations.


void pauto_pearson(
   double*       __restrict__ cov,
   double*       __restrict__ mp,
   int*          __restrict__ mpi,
   const double* __restrict__ df,
   const double* __restrict__ dg,
   const double* __restrict__ invn,
   const int tlen,
   const int offsetr,
   const int offsetc,
   const int upperbound)
{
   cov =  (double*)__builtin_assume_aligned(cov,32);
   mp =   (double*)__builtin_assume_aligned(mp,32);
   mpi =  (int*)__builtin_assume_aligned(mpi,32);
   df =   (const double*)__builtin_assume_aligned(df,32);
   dg =   (const double*)__builtin_assume_aligned(dg,32);
   invn = (const double*)__builtin_assume_aligned(invn,32);

   int rmx = (upperbound == 2*tlen) ? tlen : std::min(upperbound,tlen);
   for(int i = 0; i < rmx; i+= klen){
      int cmx = std::min(upperbound - i, tlen);
      int alignmx = cmx - cmx%klen;
      for(int j = 0; j < alignmx; j += klen){
         pauto_pearson_refkern(cov+i,mp+j,mpi+j,df+j,dg+j,invn+j,offsetr+j,offsetc+i);
      }
      if(cmx != alignmx){
         pauto_pearson_edgekern(cov+i, mp+alignmx, mpi+alignmx, df+alignmx, dg+alignmx, invn+alignmx, offsetr+alignmx, offsetc+i, cmx-alignmx);
      }
   }
}




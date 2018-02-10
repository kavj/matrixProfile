#include "../arch/v256.h"
#include<stdio.h>
#include<cmath>

typedef double dt;
// this could implement pairwise later

// It could also be optimized to do sets of 8

//template<typename dt>
void initcov(dt *a, dt *cx, dt *mu, dt *invn, dt *mp, dt *mpi, int mindiag, int len, int subLen){
   dt m1 = mu[0];
   dt s1 = invn[0];
   dt maxmp = -1;
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
      dt cor = cov * (s1 * invn[diag]);
      mp[diag] = cor;
      if(maxmp < cor){
         maxmp = cor;
      }
   }
   mpi[0] = maxmp;
   mpi[1] = maxmp;
   mpi[2] = maxmp;
   mpi[3] = maxmp;
}

typedef __m256d vt;
typedef __m256i vi;


#define stride 4

static inline void solvetile(dt*  __restrict__ cx, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dt* __restrict__ mpoff, int diag, int offset, int len) __attribute__ ((always_inline));

static inline void solvetile(dt*  __restrict__ cx, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dt* __restrict__ mpoff, int diag, int offset, int len){
   cx =  (dt*)__builtin_assume_aligned(cx,64);
   dx =  (dt*)__builtin_assume_aligned(dx,64);
   df =  (dt*)__builtin_assume_aligned(df,64);
   s  =  (dt*)__builtin_assume_aligned(s, 64);
   mp =  (dt*)__builtin_assume_aligned(mp,64);
   mpoff = (dt*)__builtin_assume_aligned(mpoff,64);
   //static unsigned long t = 0;  
   for(int i = diag; i < diag + len; i+= 32){
      vt cova, covb, covc, covd, cove, covf, covg, covh;
      strided_load_x8<vt,stride>(cova,covb,covc,covd,cove,covf,covg,covh,cx,i);
      for(int j = offset; j < offset + len; j ++){
         int k = i + j;
         vt ra = bcast<vt>(df,j);
         cova = vfma<vt>(ra,uload<vt>(dx,k),cova);
         covb = vfma<vt>(ra,uload<vt>(dx,k+stride),covb);
         covc = vfma<vt>(ra,uload<vt>(dx,k+2*stride),covc);
         covd = vfma<vt>(ra,uload<vt>(dx,k+3*stride),covd);
         cove = vfma<vt>(ra,uload<vt>(dx,k+4*stride),cove);
         covf = vfma<vt>(ra,uload<vt>(dx,k+5*stride),covf);
         covg = vfma<vt>(ra,uload<vt>(dx,k+6*stride),covg);
         covh = vfma<vt>(ra,uload<vt>(dx,k+7*stride),covh);
         
         ra = bcast<vt>(dx,j);
         cova = vfma<vt>(ra,uload<vt>(df,k),cova);
         covb = vfma<vt>(ra,uload<vt>(df,k+stride),covb);
         covc = vfma<vt>(ra,uload<vt>(df,k+2*stride),covc);
         covd = vfma<vt>(ra,uload<vt>(df,k+3*stride),covd);
         cove = vfma<vt>(ra,uload<vt>(df,k+4*stride),cove);
         covf = vfma<vt>(ra,uload<vt>(df,k+5*stride),covf);
         covg = vfma<vt>(ra,uload<vt>(df,k+6*stride),covg);
         covh = vfma<vt>(ra,uload<vt>(df,k+7*stride),covh);
      
         vt corra = cova * uload<vt>(s,k);
         vt corrb = covb * uload<vt>(s,k+stride);
         vt corrc = covc * uload<vt>(s,k+2*stride);
         vt corrd = covd * uload<vt>(s,k+3*stride);
         vt corre = cove * uload<vt>(s,k+4*stride);
         vt corrf = covf * uload<vt>(s,k+5*stride);
         vt corrg = covg * uload<vt>(s,k+6*stride);
         vt corrh = covh * uload<vt>(s,k+7*stride);
  
         ra = bcast<vt>(s,j);

         corra *= ra;
         corrb *= ra;
         corrc *= ra;
         corrd *= ra;
         corre *= ra;
         corrf *= ra;
         corrg *= ra;
         corrh *= ra;
      
         ustore<vt>(vmax(corra,uload<vt>(mp,k)),mp,k);
         ustore<vt>(vmax(corrb,uload<vt>(mp,k+stride)),mp,k+stride);
         ustore<vt>(vmax(corrc,uload<vt>(mp,k+2*stride)),mp,k+2*stride);
         ustore<vt>(vmax(corrd,uload<vt>(mp,k+3*stride)),mp,k+3*stride);
         ustore<vt>(vmax(corre,uload<vt>(mp,k+4*stride)),mp,k+4*stride);
         ustore<vt>(vmax(corrf,uload<vt>(mp,k+5*stride)),mp,k+5*stride);
         ustore<vt>(vmax(corrg,uload<vt>(mp,k+6*stride)),mp,k+6*stride);
         ustore<vt>(vmax(corrh,uload<vt>(mp,k+7*stride)),mp,k+7*stride);
         
         corra = vmax(corra,corrb);
         corrc = vmax(corrc,corrd);
         corre = vmax(corre,corrf);
         corrg = vmax(corrg,corrh);
         corra = vmax(corra,corrc);
         corre = vmax(corre,corrg);
         corra = vmax(corra,corre);
       
         astore<vt>(vmax(corra,aload<vt>(mpoff,4*j)),mpoff,4*j);
     //    t+=64;
      }
      strided_store_x8<vt,stride>(cova,covb,covc,covd,cove,covf,covg,covh,cx,i);
   }
   //printf("%lu\n",t);
}

void solvempref(dt* __restrict__ a, dt*  __restrict__ cx, dt* __restrict__ mu, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dt* __restrict__ mpi, int len, int diagmin, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dt*)__builtin_assume_aligned(mpi,64);
   const int tilelen = 512;
   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
   int aldiagmax = len - sublen - tilelen -
   for(int i = diagmin; i < len-sublen-tilelen; i+=tilelen){
      int jmax = len-i-tilelen;
      for(int j = 0; j < jmax; j+= tilelen){
         solvetile(cx,df,dx,s,mp,mpi,i,j,tilelen);
      }
   }
}



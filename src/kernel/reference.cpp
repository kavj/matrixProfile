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
      dt cor = cov * (s1 * invn[diag]);
      mp[diag] = cor;
      mpi[diag] = 0;
      if(maxmp < cor){
         maxmp = cor;
         maxmpi = diag;
      }
   }
   mp[0] = maxmp;
   mpi[0] = maxmpi;
}


/*
static inline void solvetile(dt*  __restrict__ cx, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int diag, int offset, int len) __attribute__ ((always_inline));

static inline void solvetile(dt*  __restrict__ cx, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int diag, int offset, int len){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);
   for(int i = diag; i < diag + len; i++){
      dt cov = cx[i];
      for(int j = offset; j < offset + len; j++){
         int k = i+j;
         cov = fma(df[j],dx[i+j],cov); 
         cov = fma(df[i+j],dx[j],cov);
         dt cor = cov*s[j];
         cor *= s[k];
         dt mpa = mp[j];
         dt mpb = mp[k];
         if(mpa < cor){
            mp[j] = cor;
            mpi[j] = k;
         }
         if(mpb < cor){
            mp[k] = cor;
            mpi[k] = j;
         }
      }
   }
   
}
*/


static inline void solvetile(dt* __restrict__cx, const dt* __restrict__ df, const dt* __restrict__dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len) __attribute__ ((always_inline));
//solvempref(T,cx,mu,dF,dX,sigmaInv,buffer,bufferI,n,m,m);

static inline void solvetile(dt* __restrict__ cx, const dt* __restrict__ df, const dt* __restrict__ dx, const dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int dm, int om, int len){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);

   for(int diag = dm; diag < dm+8; diag+=8){
      dt c1 = cx[diag];
      dt c2 = cx[diag+1];
      dt c3 = cx[diag+2];
      dt c4 = cx[diag+3];
      dt c5 = cx[diag+4];
      dt c6 = cx[diag+5];
      dt c7 = cx[diag+6];
      dt c8 = cx[diag+7];
      for(int offset = om; offset < om+len; offset++){
         int k = offset+diag;
         dt ea = dx[offset];
         c1 = fma(ea,df[k],c1);
         c2 = fma(ea,df[k+1],c2);
         c3 = fma(ea,df[k+2],c3);
         c4 = fma(ea,df[k+3],c4);
         c5 = fma(ea,df[k+4],c5);
         c6 = fma(ea,df[k+5],c6);
         c7 = fma(ea,df[k+6],c7);
         c8 = fma(ea,df[k+7],c8);

         ea = df[offset];
         c1 = fma(ea,dx[k],c1);
         c2 = fma(ea,dx[k+1],c2);
         c3 = fma(ea,dx[k+2],c3);
         c4 = fma(ea,dx[k+3],c4);
         c5 = fma(ea,dx[k+4],c5);
         c6 = fma(ea,dx[k+5],c6);
         c7 = fma(ea,dx[k+6],c7);
         c8 = fma(ea,dx[k+7],c8);
         
         dt cr1 = c1 * s[k];
         dt cr2 = c2 * s[k+1];
         dt cr3 = c3 * s[k+2];
         dt cr4 = c4 * s[k+3];
         dt cr5 = c5 * s[k+4];
         dt cr6 = c6 * s[k+5];
         dt cr7 = c7 * s[k+6];
         dt cr8 = c8 * s[k+7];

         ea = s[offset];

         cr1 *= ea;
         cr2 *= ea;
         cr3 *= ea;
         cr4 *= ea;
         cr5 *= ea;
         cr6 *= ea;
         cr7 *= ea; 
         cr8 *= ea;
 
         if(mp[k] < cr1){
            mp[k] = cr1;
            mpi[k] = offset;
         }
         if(mp[k+1] < cr2){
            mp[k+1] = cr2;
            mpi[k+1] = offset;
         }
         if(mp[k+2] < cr3){
            mp[k+3] = cr3;
            mpi[k+3] = offset;
         } 
         if(mp[k+4] < cr4){
            mp[k+4] = cr4;
            mpi[k+4] = offset;
         }
         if(mp[k+5] < cr5){
            mp[k+5] = cr5;
            mpi[k+5] = offset;
         }
         if(mp[k+6] < cr6){
            mp[k+6] = cr6;
            mpi[k+6] = offset;
         }
         if(mp[k+7] < cr7){
            mp[k+7] = cr7;
            mpi[k+7] = offset;
         }
         if(mp[k+8] < cr8){
            mp[k+8] = cr8;
            mpi[k+8] = offset;
         }

         int sel1 = k; 
         int sel2 = k+2;
         int sel3 = k+4;
         int sel4 = k+6;
          
         if(cr1 < cr2){
            cr1 = cr2;
            sel1 = k+1;
         }
         if(cr3 < cr4){
            cr3 = cr4;
            sel2 = k+2;
         }
         if(cr5 < cr6){
            cr5 = cr6;
            sel3 = k+5;
         }
         if(cr7 < cr8){
            cr7 = cr8;
            sel4 = k+7;
         }
         if(cr1 < cr3){
            cr1 = cr3; 
            sel1 = sel2;
         }
         if(cr5 < cr7){
            cr5 = cr7;
            sel3 = sel4;
         }
         if(cr1 < cr5){
            cr1 = cr5;
            sel1 = sel3;
         }
         if(cr1 < mp[offset]){
            mp[offset] = cr1;
            mpi[offset] = sel1;
         }

      }
      cx[diag] = c1;
      cx[diag+1] = c2;
      cx[diag+2] = c3;
      cx[diag+3] = c4;
      cx[diag+4] = c5;
      cx[diag+5] = c6;
      cx[diag+6] = c7;
      cx[diag+7] = c8;
 
   }
}


/*
void solvempref(dt* __restrict__ a, dt*  __restrict__ cx, dt* __restrict__ mu, dt* __restrict__ df, dt* __restrict__ dx, dt* __restrict__ s, dt* __restrict__ mp, dti* __restrict__ mpi, int len, int diagmin, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);
   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
   const int tilewid = 512;
   for(int i = diagmin; i < len-sublen-tilewid; i+= 8){
      int jmax = len - sublen - i + 1 - tilewid;
      //for(int j = 1; j < jmax; j+= tilewid){
         solvetile(cx,df,dx,s,mp,mpi,i,1,jmax);
   }
}
*/


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
      for(int j = 1; j < jmax; j++){
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


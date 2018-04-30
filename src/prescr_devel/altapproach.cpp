#include "../arch/avx256.h"
#include "../utils/reg.h"
using namespace avx256_t;

/*
void bleh2(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
      
   for(int j = 0; j < 64; j++){
      cov[j] += dx[j]*df[j+offsetc];
      cov[j] += df[j]*dx[j+offsetc];
      if(mp[j] < cov[j]*s[j]*s[j+offsetc]){
         mp[j] = cov[j]*s[j]*s[j+offsetc];
      }
      if(mp[j+offsetc] < cov[j]*s[j]*s[j+offsetc]){
         mp[j+offsetc] = cov[j]*s[j]*s[j+offsetc];
      }
   }
}
*/

/*
void  bleh(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);
   for(int j = 0; j < 64; j++){
      cov[j] += dx[j]*df[j+offsetc];
      cov[j] += df[j]*dx[j+offsetc];
      if(mp[j] < cov[j]*s[j]*s[j+offsetc]){
         mp[j] = cov[j]*s[j]*s[j+offsetc];
         mpi[j] = offsetr+j;
      }
      if(mp[j+offsetc] < cov[j]*s[j]*s[j+offsetc]){
         mp[j+offsetc] = cov[j]*s[j]*s[j+offsetc];
         mpi[j+offsetc] = j;
      }
   }
}*/


#define tsz 64
#define unroll 8
void  bleh4(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);

   for(int i = 0; i < tsz; i++){
      for(int j = 0; j < tsz; j+=32){
         block<__m256d> cov_r;
         for(int k = 0; k < unroll; k++){
            cov_r(k) = aload(cov,4*(j+k));
         }
         __m256d q = brdcst(dx,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(df,i+j+4*k),cov_r(k));
         }
         q = brdcst(df,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) = mul_add(q,uload(dx,i+j+4*k),cov_r(k));
            astore(cov_r(k),cov,4*(j+k));
         }
         q = brdcst(s,i);
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= q;
         }
         for(int k = 0; k < unroll; k++){
            cov_r(k) *= uload(s,i+j+4*k);
         }
         block<__m256i> mask;
         for(int k = 0; k < unroll; k++){
            mask(k) = cov_r(k) > uload(mp,i+j+4*k);
         }
         for(int k = 0; k < unroll; k+=2){
            __m256i r = brdcst(i);
            int l = testnz(mask(k));
            int m = testnz(mask(k+1));
            if(l){
               if(m){
                  maskstore(r,mask(k),mpi+i+j+4*k);
                  maskstore(r,mask(k+1),mpi+i+j+4*k+4);
                  maskstore(cov_r(k),mask(k),mp+i+j+4*k);
                  maskstore(cov_r(k+1),mask(k+1),mp+i+j+4*k+4);
               }
               else{
                  maskstore(cov_r(k),mask(k),mp+i+j+4*k);
                  maskstore(r,mask(k),mpi+i+j+4*k+4);
               }
            }
            else if(m){
               maskstore(cov_r(k),mask(k),mp+i+j+4*k);
               maskstore(r,mask(k),mpi+i+j+4*k+4);
            }
         }
         for(int k = 0; k < unroll/2; k++){
            mask(k) = cov_r(2*k) > cov_r(2*k+1);
         }
         __m256i base = set(0,1,2,3);
         for(int k = 0; k < unroll/2; k++){
            if(testnz(mask(2*k))){
               mask(k+unroll/2) = blend(
            }
         }
         for(int k = 0; k < unroll/4; k++){
         
         }
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i] < cov[j]*s[i]*s[i+j+offsetc]){
            mp[i] = cov[j]*s[i]*s[i+j+offsetc];
            mpi[i] = i+offsetc+j;
         }
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i+j+offsetc] < cov[j]*s[i]*s[i+j+offsetc]){
            mp[i+j+offsetc] = cov[j]*s[i]*s[i+j+offsetc];
            mpi[i+j+offsetc] = i+offsetc+j;
         }
      }*/
      }
   }
}





#define tsz 64 
void  bleh(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);

   for(int i = 0; i < tsz; i++){
      for(int j = 0; j < tsz; j++){
         cov[j] += dx[i]*df[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         cov[j] += df[i]*dx[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i] < cov[j]*s[i]*s[i+j+offsetc]){
            mp[i] = cov[j]*s[i]*s[i+j+offsetc];
            mpi[i] = i+offsetc+j;
         }
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i+j+offsetc] < cov[j]*s[i]*s[i+j+offsetc]){
            mp[i+j+offsetc] = cov[j]*s[i]*s[i+j+offsetc];
            mpi[i+j+offsetc] = i+offsetc+j;
         }
      }
   }
}


/*
void  bleh_avx(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);

   int unroll = 32;
   
   for(int i = 0; i < tsz; i++){
      block<__m256d> cov_r;
      for(int j = 0; j < tsz; j+=unroll){
         cov_r = aload(cov,4*j);
      }
      __m256d q = brdcst(dx,;
      for(int j = 0; j < tsz; j+=unroll){
         cov_r = 
      }
      block<__m256d> cov_r;
      __m256d dx_r = brdcst(dx,i+j);
      for(int j = 0; j < unroll; i++){
         cov_r(i) = mul_add(dx_r,uload(df,i+j+
      }  
      for(int j = 0; j < tsz; j++){
         cov[j] += dx[i+j]*df[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         cov[j] += df[i+j]*dx[i+j+offsetc];
      }
      for(int j = 0; j < tsz; j++){
         if(mp[i+j] < cov[j]*s[i+j]*s[i+j+offsetc]){
            mp[i+j] = cov[j]*s[i+j]*s[i+j+offsetc];
            mpi[i+j] = offsetr+j;
         }
         if(mp[i+j+offsetc] < cov[j]*s[j]*s[j+offsetc]){
            mp[i+j+offsetc] = cov[j]*s[j]*s[j+offsetc];
            mpi[i+j+offsetc] = j;
         }
      }
   }
}
*/






/*
void bleh(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);
      
   for(int j = 0; j < 64; j++){
      cov[j] += dx[j]*df[j+offsetc];
      cov[j] += df[j]*dx[j+offsetc];
      if(mp[j] < cov[j]*s[j]*s[j+offsetc]){
         mp[j] = cov[j]*s[j]*s[j+offsetc];
         mpi[j] = offsetr+j;
      }
      if(mp[j+offsetc] < cov[j]*s[j]*s[j+offsetc]){
         mp[j+offsetc] = cov[j]*s[j]*s[j+offsetc];
         mpi[j+offsetc] = j;
      }
   }
}
*/


/*void bleh(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);
      
   double corr[64];
   for(int j = 0; j < 64; j++){
      cov[j] += dx[j]*df[j+offsetc];
   }
   for(int j = 0; j < 64; j++){
      cov[j] += df[j]*dx[j+offsetc];
      corr[j] = cov[j]*s[j]*s[j+offsetc];
   }
   for(int j = 0; j < 64; j++){
      if(corr[j] < mp[j]){
         mp[j] = corr[j];
         mpi[j] = offsetr+j;
      }
      if(corr[j] < mp[j+offsetc]){
         mp[j+offsetc] = corr[j];
         mpi[j+offsetc] = j+offsetc;
      }
   }
}*/


/*
void bleh2(double* __restrict__ cov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ s, double* __restrict__ mp, int* __restrict__ mpi, int offsetr, int offsetc, int offsetmp){
   cov = (double*)__builtin_assume_aligned(cov,32);
   df = (const double*)__builtin_assume_aligned(df,32);
   dx = (const double*)__builtin_assume_aligned(dx,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi = (int*)__builtin_assume_aligned(mpi,32);
      
   block<__m256d> r;
   __m256d q = brdcst(dx,offsetc);
   for(int j = 0; j < 8; j++){
      r(j) = mul_add(q,aload(df,0),aload(cov,4*j));
   }
   q = brdcst(df,offsetc);
   for(int j = 0; j < 8; j++){
      r(j) = mul_add(q,aload(dx,0),r(j));
   }
   q = brdcst(s,offsetc);
   for(int j = 0; j < 8; j++){
      astore(r(j),cov,4*j);
      r(j) *= q;
   }
   for(int j = 0; j < 8; j++){
      r(j) *= uload(s,j);
   }
   for(int j = 0; j < 8; j++){
      r(j) = vmax(r(j),aload(mp,4*j));
      astore(r(j),mp,4*j);
   }
   
}

*/



/*
void bleh_long(const double* __restrict__ a, const double* __restrict__ b, const double* __restrict__ e, double* __restrict__ c, long* __restrict__ d, int len, int k){
   a = (const double*)__builtin_assume_aligned(a,64);
   b = (const double*)__builtin_assume_aligned(b,64);
   c = (double*)__builtin_assume_aligned(c,64);
   d = (long*)__builtin_assume_aligned(d,64);
   e = (const double*)__builtin_assume_aligned(e,64);
   for(int j = 0; j < 64; j++){
   for(int i = 0; i < 64; i++){
      c[i] += a[i+j]*b[j+k];
   }
   }

*/
  /* a = (double*)__builtin_assume_aligned(a,32);
   b = (double*)__builtin_assume_aligned(b,32);
   c = (double*)__builtin_assume_aligned(c,32);
   d = (int*)__builtin_assume_aligned(d,32);


   for(int i = 0; i < 32; i++){
      for(int j = 0; j < 16; j++){
  //       c[i] = a[i+j] + b[i+j]*e[i];    

         if(a[i] < b[i]){
            c[i] = b[i];
            d[i] = j+k;
         }
      }
  
   }
}
*/
/*
void extrap(double* __restrict__ cov, double* __restrict__ mp, long* __restrict mpi, const double* __restrict__ s, const double* __restrict__ df, const double* __restrict__ dg, const double* __restrict__ dh, const double* __restrict dx,int offset){
   cov = (double*) __builtin_assume_aligned(cov,32);
   s = (const double*)__builtin_assume_aligned(s,32);
   df = (const double*)__builtin_assume_aligned(df,32); 
   dx = (const double*)__builtin_assume_aligned(dx,32);
   dg = (const double*)__builtin_assume_aligned(dg,32);
   dh = (const double*)__builtin_assume_aligned(dh,32);
   mp = (double*)__builtin_assume_aligned(mp,32);
   mpi= (long*)__builtin_assume_aligned(mpi,32);
    
   for(int l = 0; l < 32; l++){ 
   for(int i = l; i < l+64; i++){
      for(int j = 0; j < 8; j++){
         cov[i+j] += df[i+j];
      }
      for(int j = 0; j < 16; j++){
         cov[i+j] += dg[i+j]*dh[j];
      }
      
   }
   }
}


*/


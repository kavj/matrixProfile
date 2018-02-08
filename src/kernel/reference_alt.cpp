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


//solvempref(T,cx,mu,dF,dX,sigmaInv,buffer,bufferI,n,m,m);




void solvempref(dt *a, dt *cx, dt *mu, dt *df, dt *dx, dt *s, dt *mp, dti *mpi, int len, int diagmin, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   dx = (dt*)__builtin_assume_aligned(dx,64);
   df = (dt*)__builtin_assume_aligned(df,64);
   s  = (dt*)__builtin_assume_aligned(s,64);
   mp =(dt*)__builtin_assume_aligned(mp,64);
   mpi = (dti*)__builtin_assume_aligned(mpi,64);

   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
/*   for(int i = 0; i < 10000; i+=100){
      printf("%lf %d \n",cx[i],i);
   }
*/
   for(int i = diagmin; i < len-sublen; i++){
      int jmax = len - sublen - i + 1;
      dt cov = cx[i];
      for(int j = 1; j < jmax; j++){
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
      cx[i] = cov; 
   }
}


#include<stdio.h>
#include<cmath>

typedef double dt;
typedef long dti;
// this could implement pairwise later

// It could also be optimized to do sets of 8

//template<typename dt>
void initcov(dt *a, dt *mu, dt *invn, dt *c, dt *mp, dti *mpi, int mindiag, int len, int subLen){
   dt m1 = mu[0];
   dt s1 = invn[0];
   dt maxmp = -1;
   dti maxmpi = -1;
   int upperlim = len - subLen;
   for(int diag = mindiag; diag < upperlim; diag++){
      dt cov = 0;
      dt m2 = mu[diag];
      for(int offset = 0; offset < subLen; offset++){
         dt c1 = a[offset] - m1;
         dt c2 = a[diag+offset] - m2;
         cov = fma(c1,c2,cov);
      }
      c[diag] = cov;
      dt cor = cov * (s1 * invn[diag]);
      mp[diag] = cor;
      mpi[diag] = 1;
      printf("cov: %lf, cor %lf diag %d\n",cov,cor,diag);
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
   initcov(a,cx,mu,s,mp,mpi,diagmin,len,sublen);
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
         if(j < diagmin){
            printf("%lf %lf %d\n",mpa,cor,j);
         }
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


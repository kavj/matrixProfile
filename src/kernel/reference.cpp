#include<cmath>

typedef double dt;
typedef long dti;
// this could implement pairwise later

// It could also be optimized to do sets of 8

//template<typename dt>
void init(dt *a, dt *mu, dt *invn, dt *c, dt *mp, dti *mpi, int mindiag, int len, int subLen){
   dt m1 = mu[0];
   dt s1 = invn[0];
   int upperlim = len - subLen;
   for(int diag = mindiag; diag < upperlim; diag++){
      dt cov = 0;
      dt m2 = mu[diag];
      for(int offset = 0; offset < subLen; offset++){
         double c1 = a[offset] - m1;
         double c2 = a[diag+offset] - m2;
         cov = fma(c1,c2,cov);
      }
      c[diag] = cov * (s1 * invn[diag]);
   }
}




void solvempref(dt *cx, dt *mu, dt *df, dt *dx, dt *s, dt *mp, dti *mpi, int len, int diagmin, int sublen){
   for(int i = diagmin; i < diagmax; i++){
      int jmax = len - sublen - i;
      for(int j = 1; j < jmax; j++){
         int k = i+j;
         dt cov = fma(df[j],dg[i+j],c[i]); 
         cor = fma(df[i+j],dg[j],cov);
         c[i] = cov;
         cor *= s[j];
         cor *= s[k];
         mpa = mp[j];
         mpb = mp[k];
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




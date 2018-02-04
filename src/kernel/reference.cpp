

template<typename dt>
void init(dt *a, dt *mu, dt *invn, dt *c, int mindiag, int len, int subLen){
   dt m2 = mu[0];
   dt s1 = invn[0];
   int upperlim = len - subLen; 
   for(int diag = mindiag; diag < upperlim; diag++){
      dt c1 = 0; 
      dt m1 = mu[diag];
      for(int offset = 0; offset < subLen; offset++){
         c1 += (a[diag+offset] - m1) * (a[offset] - m2);
      }
      c[diag] = c1 * (s1 * invn[diag]);
   }
}





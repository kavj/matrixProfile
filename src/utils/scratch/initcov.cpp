#include "../arch/avx256.h"
//#include<immintrin.h>
#include <cstdlib>
#include "reg.h"

using namespace avx256_t;
typedef double dt;


void initcov(const dt* __restrict__ a, dt* __restrict__ cx, const dt* __restrict__ mu, int mindiag, int len, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   a  = (dt*)__builtin_assume_aligned(a,64);
   mu = (dt*)__builtin_assume_aligned(mu,64);
 
   dt m1 = mu[0];
   int upperlim = len - sublen + 1;
   int allim = upperlim - upperlim%32;
   dt* q = static_cast<double*>(malloc(sublen*sizeof(dt)));
   for(int i = 0; i < sublen; i++){
      q[i] = a[i] - m1;
   }
   for(int diag = mindiag; diag < allim; diag+=32){
      block<__m256d> cov;
      for(int offset = 0; offset < sublen; offset++){
         __m256d qi = brdcst(q,offset);
         for(int subd = 0; subd < 8; subd++){
            cov(subd+8) = uload(a,diag+offset+subd) - uload(mu,diag+subd);
         }
         for(int subd = 0; subd < 8; subd++){
            cov(subd) = mul_add(qi,cov(subd+8),cov(subd));
         }
      }
      for(int subd = 0; subd < 8; subd++){
         astore(cov(subd),cx,diag+4*subd);
      }
   }
   for(int diag = allim; diag < upperlim; diag++){
      dt covs = 0;
      for(int offset = 0; offset < sublen; offset++){
         covs += (a[diag+offset] - mu[diag])*q[offset];
      }
   }
   free(q);
}


void initcov_brdcst(const dt* __restrict__ a, dt* __restrict__ cx, const dt* __restrict__ mu, int mindiag, int len, int sublen){
   cx = (dt*)__builtin_assume_aligned(cx,64);
   a  = (dt*)__builtin_assume_aligned(a,64);
   mu = (dt*)__builtin_assume_aligned(mu,64);
 
   dt m1 = mu[0];
   int upperlim = len - sublen + 1;
   int allim = upperlim - upperlim%32;
   dt* q = static_cast<double*>(malloc(sublen*sizeof(dt)));
   for(int i = 0; i < sublen; i++){
      q[i] = a[i] - m1;
   }
   for(int diag = mindiag; diag < allim; diag+=32){
      block<__m256d> cov;
      for(int offset = 0; offset < sublen; offset++){
         __m256d qi = brdcst(q,offset);
         for(int subd = 0; subd < 8; subd++){
            cov(subd+8) = uload(a,diag+offset+subd) - uload(mu,diag+subd);
         }
         for(int subd = 0; subd < 8; subd++){
            cov(subd) = mul_add(qi,cov(subd+8),cov(subd));
         }
      }
      for(int subd = 0; subd < 8; subd++){
         astore(cov(subd),cx,diag+4*subd);
      }
   }
   for(int diag = allim; diag < upperlim; diag++){
      dt covs = 0;
      for(int offset = 0; offset < sublen; offset++){
         covs += (a[diag+offset] - mu[diag])*q[offset];
      }
   }
   free(q);
}



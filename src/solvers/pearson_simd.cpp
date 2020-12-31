#include<immintrin.h>


/*

Note to colleagues, who may be reading this.
This has no knowledge of how to handle nans. 
It is optimized for finite data. Intel implements max

as  max(a,b) = a > b ? a : b

which makes these sensitive to reordering.

There are suggested ways to hack around this, but they're very inadvisable.

The arcane names refer to row/column indexing with respect to an implicitly defined correlation matrix.

I used this notation in matlab previously, even though it's not ideal.

The fwd designation indicates that it's used to add the product of a pair of elements to a running sum.
The bwd designation indicates that it's used to remove the product of a trailing pair.

This formula tends to blow up less on singularities (streams of constants), but I still suggest preconditioning.

I included a basic preconditioner with the Matlab front end.


Since this is specialized for AVX2, it computes sections of 32 "diagonals" at a time.
"iters" is the number of "rows" traversed.

Since this problem has skewed boundary conditions, "iters" must take that into consideration.
That is to say, boundary conditions here are set by the most restrictive component.

mp and invn take both row and column views. 


*/


void pearson_inner(
   double*       __restrict cv,
   double*       __restrict mp,
   const double* __restrict dr_bwd,
   const double* __restrict dc_bwd,
   const double* __restrict dr_fwd,
   const double* __restrict dc_fwd,
   const double* __restrict invn,
   const int colOf,
   const int iters)
{

   __m256d c0 = _mm256_load_pd(cv);
   __m256d c1 = _mm256_load_pd(cv + 4);
   __m256d c2 = _mm256_load_pd(cv + 8);
   __m256d c3 = _mm256_load_pd(cv + 12);
   __m256d c4 = _mm256_load_pd(cv + 16);
   __m256d c5 = _mm256_load_pd(cv + 20);
   __m256d c6 = _mm256_load_pd(cv + 24);
   __m256d c7 = _mm256_load_pd(cv + 28);
   for(int i = 0; i < iters; i++){
      int j = i - 1;
      if(i > 0){
         __m256d dr_bwd = _mm256_broadcast_sd(dr_bwd + j);

         // cv - dr_bwd[r] * dc_bwd[col]
         c0 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j), c0);
         c1 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 4), c1);
         c2 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 8), c2);
         c3 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 12), c3);
         c4 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 16), c4);
         c5 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 20), c5);
         c6 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 24), c6);
         c7 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + j + 28), c7);

         __m256d dr_fwd_ = _mm256_broadcast_sd(dr_fwd + j);

         // cv + dr_fwd[r] * dc_fwd[col]
         c0 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j), c0);
         c1 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 4), c1);
         c2 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 8), c2);
         c3 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 12), c3);
         c4 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 16), c4);
         c5 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 20), c5);
         c6 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 24), c6);
         c7 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + j + 28), c7);
      }
      int k = i + colOf;
      __m256d invn_ = _mm256_broadcast_sd(invn + i);
      __m256d f0 = _mm256_mul_pd(c0, _mm256_loadu_pd(invn + k));
      __m256d f1 = _mm256_mul_pd(c1, _mm256_loadu_pd(invn + k + 4));
      __m256d f2 = _mm256_mul_pd(c2, _mm256_loadu_pd(invn + k + 8));
      __m256d f3 = _mm256_mul_pd(c3, _mm256_loadu_pd(invn + k + 12));
      __m256d f4 = _mm256_mul_pd(c4, _mm256_loadu_pd(invn + k + 16));
      __m256d f5 = _mm256_mul_pd(c5, _mm256_loadu_pd(invn + k + 20));
      __m256d f6 = _mm256_mul_pd(c6, _mm256_loadu_pd(invn + k + 24));
      __m256d f7 = _mm256_mul_pd(c7, _mm256_loadu_pd(invn + k + 28));

      __m256d g0 = _mm256_mul_pd(f0, invn_);
      __m256d g1 = _mm256_mul_pd(f1, invn_);
      __m256d g2 = _mm256_mul_pd(f2, invn_);
      __m256d g3 = _mm256_mul_pd(f3, invn_);
      __m256d g4 = _mm256_mul_pd(f4, invn_);
      __m256d g5 = _mm256_mul_pd(f5, invn_);
      __m256d g6 = _mm256_mul_pd(f6, invn_);
      __m256d g7 = _mm256_mul_pd(f7, invn_);

      __m256d h0 = _mm256_max_pd(g0, _mm256_loadu_pd(mp + k));
      __m256d h1 = _mm256_max_pd(g1, _mm256_loadu_pd(mp + k + 4));
      __m256d h2 = _mm256_max_pd(g2, _mm256_loadu_pd(mp + k + 8));
      __m256d h3 = _mm256_max_pd(g3, _mm256_loadu_pd(mp + k + 12));
      __m256d h4 = _mm256_max_pd(g4, _mm256_loadu_pd(mp + k + 16));
      __m256d h5 = _mm256_max_pd(g5, _mm256_loadu_pd(mp + k + 20));
      __m256d h6 = _mm256_max_pd(g6, _mm256_loadu_pd(mp + k + 24));
      __m256d h7 = _mm256_max_pd(g7, _mm256_loadu_pd(mp + k + 28));

      _mm256_storeu_pd(mp + k, h0);
      _mm256_storeu_pd(mp + k + 4, h1);
      _mm256_storeu_pd(mp + k + 8, h2);
      _mm256_storeu_pd(mp + k + 12, h3);
      _mm256_storeu_pd(mp + k + 16, h4);
      _mm256_storeu_pd(mp + k + 20, h5);
      _mm256_storeu_pd(mp + k + 24, h6);
      _mm256_storeu_pd(mp + k + 28, h7);

      __m256d j0 = _mm256_max_pd(g0, g1);
      __m256d j1 = _mm256_max_pd(g2, g3);
      __m256d j2 = _mm256_max_pd(g4, g5);
      __m256d j3 = _mm256_max_pd(g6, g7);

      __m256d L0 = _mm256_max_pd(j0, j1);
      __m256d L1 = _mm256_max_pd(j2, j3);
   
      __m256d L2 = _mm256_max_pd(L0, L1);

      double q[4];
      _mm256_store_pd(&q[0], L2);
      double q0 = q[0] > q[1] ? q[0] : q[1];
      double q1 = q[2] > q[3] ? q[2] : q[3];
      double q2 = q0 > q1 ? q0 : q1;
      mp[i] = q2 > mp[i] ? q2 : mp[i];
   }
   /*
     These stores may be dead if you're not doing any further
     processing on this block. They are however useful if you 
     need to debug divergence issues, by comparing to a reference
     version of what you would expect for the last value.

   */
   _mm256_store_pd(cv, c0);
   _mm256_store_pd(cv + 4, c1);
   _mm256_store_pd(cv + 8, c2);
   _mm256_store_pd(cv + 12, c3);
   _mm256_store_pd(cv + 16, c4);
   _mm256_store_pd(cv + 20, c5);
   _mm256_store_pd(cv + 24, c6);
   _mm256_store_pd(cv + 28, c7);
}


/*

  This is only used for constrained cases. 
  Otherwise we use the simd version.

  This assumes that at this level, cov[d] is active as long as
  d < clim - r

*/

static inline void pearson_edge(
   double*       __restrict cov, 
   double*       __restrict mp,  
   const double* __restrict dr_bwd,
   const double* __restrict dc_bwd,
   const double* __restrict dr_fwd,
   const double* __restrict dc_fwd,
   const double* __restrict invn,
   int clim)
{
   for(int r = 0; r < clim; r++){
      for(int d = 0; d < clim - r; d++){
         int col = r + d;
         if(r > 0){
            cov[d] -= dr_bwd[r-1] * dc_bwd[col-1];
            cov[d] += dr_fwd[r-1] * dc_fwd[col-1];
         }
         double cr = cov[d] * invn[r] * invn[col];
         if(mp[r] < cr){
            mp[r] = cr;
         }
         if(mp[col] < cr){
            mp[col] = cr;
         }
      }
   }
}


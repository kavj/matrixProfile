#include <immintrin.h>
#include <algorithm>

/*

I am leaving a few notes for colleagues, who may be reading this.


First, regarding nans and subsequences containing only constants.

The methods at this level have no knowledge of these things.
They are optimized for finite data, which has a couple edge cases.


Since not everyone is aware of this, I'll point out the source of divergence in cases where intermediate arithmetic produces nan values
(typically via overflow followed by multiplication)

The most common assembly implementation of max (which is not compliant with IEEE2018 rules)

is  

max(a,b) = a > b ? a : b


This means:

if a or b are +- nan, then 

a > b == false 

and so

max(a,b) = b 

Now if you have a cascaded reduction

max(<a0,b0>.....<ak,bk>)

it's possible for the largest value to be in a comparison of the form

max(largest, nan) --> nan

with the next comparison being

max(nan, some_other_value)

in which case, some_other_value is silently propagated

IEEE2018 prefers to propagate nan. Matlab prefers to propagate a finite value, which is error prone.


Second note...

I have made considerable efforts to improve the stability of lookahead methods wherever possible. 
I swapped out the method used for difference equations.
I provided some preprocessing to slightly perturb constant regions.

The last point will be to explicitly detect that divergence has occurred. I used this in writing unit tests previously
and suggested it to others before that point, but I'm now embedding it in the calculations.

Since it's cheap enough to naively compute a linear number of dot products, we can do that to explicitly check the final
output of the variable "cv" against reference values, computed naively.

If any of these diverge, then we issue a warning and possibly output the corresponding "diagonal" indices.

Last note on implementation...

This is specialized for AVX2. If I include index calculations and want to maintain a reasonably fast implementation,
I have to use masked writes to minimize register spill. Masked writes do not have good throughput on AMD hardware. Seeing as I
prefer not to deal with this, I have dropped index calculations.

This implementation is meant to replace the use of "anytime" methods, where you don't get an exact index anyway.
If you just want the neighbors of the low values, they can be cheaply computed in isolation.


Reference
https://www.agner.org/optimize/instruction_tables.pdf


*/



static inline void pearson_inner(
    double* __restrict cv,
    double* __restrict mpr,
    double* __restrict mpc,
    const double* __restrict dr_bwd,
    const double* __restrict dc_bwd,
    const double* __restrict dr_fwd,
    const double* __restrict dc_fwd,
    const double* __restrict invnr,
    const double* __restrict invnc,
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
   for (int i = 0; i < iters; i++)
   {
      if (i > 0)
      {
         // The difference equation indices fold a minus 1 calculation into the constant offset. 
         // It's attributable to a maximum of subsequence count - 1 steps in total for subsequence count windows
         __m256d dr_bwd_ = _mm256_broadcast_sd(dr_bwd + i - 1);

         // cv - dr_bwd[r] * dc_bwd[col]
         c0 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i - 1), c0);
         c1 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 3), c1);
         c2 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 7), c2);
         c3 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 11), c3);
         c4 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 15), c4);
         c5 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 19), c5);
         c6 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 23), c6);
         c7 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 27), c7);

         __m256d dr_fwd_ = _mm256_broadcast_sd(dr_fwd + i - 1);

         // cv + dr_fwd[row] * dc_fwd[col]
         c0 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i - 1), c0);
         c1 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 3), c1);
         c2 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 7), c2);
         c3 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 11), c3);
         c4 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 15), c4);
         c5 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 19), c5);
         c6 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 23), c6);
         c7 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 27), c7);
      }
      // invn and mp don't need to fold a -1, so their constants are 1 step higher
      __m256d invn_ = _mm256_broadcast_sd(invnr + i);
      __m256d f0 = _mm256_mul_pd(c0, _mm256_loadu_pd(invnc + i));
      __m256d f1 = _mm256_mul_pd(c1, _mm256_loadu_pd(invnc + i + 4));
      __m256d f2 = _mm256_mul_pd(c2, _mm256_loadu_pd(invnc+ i + 8));
      __m256d f3 = _mm256_mul_pd(c3, _mm256_loadu_pd(invnc+ i + 12));
      __m256d f4 = _mm256_mul_pd(c4, _mm256_loadu_pd(invnc + i + 16));
      __m256d f5 = _mm256_mul_pd(c5, _mm256_loadu_pd(invnc + i + 20));
      __m256d f6 = _mm256_mul_pd(c6, _mm256_loadu_pd(invnc + i + 24));
      __m256d f7 = _mm256_mul_pd(c7, _mm256_loadu_pd(invnc + i + 28));

      __m256d g0 = _mm256_mul_pd(f0, invn_);
      __m256d g1 = _mm256_mul_pd(f1, invn_);
      __m256d g2 = _mm256_mul_pd(f2, invn_);
      __m256d g3 = _mm256_mul_pd(f3, invn_);
      __m256d g4 = _mm256_mul_pd(f4, invn_);
      __m256d g5 = _mm256_mul_pd(f5, invn_);
      __m256d g6 = _mm256_mul_pd(f6, invn_);
      __m256d g7 = _mm256_mul_pd(f7, invn_);

      __m256d h0 = _mm256_max_pd(g0, _mm256_loadu_pd(mpc + i));
      __m256d h1 = _mm256_max_pd(g1, _mm256_loadu_pd(mpc + i + 4));
      __m256d h2 = _mm256_max_pd(g2, _mm256_loadu_pd(mpc + i + 8));
      __m256d h3 = _mm256_max_pd(g3, _mm256_loadu_pd(mpc + i + 12));
      __m256d h4 = _mm256_max_pd(g4, _mm256_loadu_pd(mpc + i + 16));
      __m256d h5 = _mm256_max_pd(g5, _mm256_loadu_pd(mpc + i + 20));
      __m256d h6 = _mm256_max_pd(g6, _mm256_loadu_pd(mpc + i + 24));
      __m256d h7 = _mm256_max_pd(g7, _mm256_loadu_pd(mpc + i + 28));

      _mm256_storeu_pd(mpc, h0);
      _mm256_storeu_pd(mpc + 4, h1);
      _mm256_storeu_pd(mpc + 8, h2);
      _mm256_storeu_pd(mpc + 12, h3);
      _mm256_storeu_pd(mpc + 16, h4);
      _mm256_storeu_pd(mpc + 20, h5);
      _mm256_storeu_pd(mpc + 24, h6);
      _mm256_storeu_pd(mpc + 28, h7);

      __m256d j0 = _mm256_max_pd(g0, g1);
      __m256d j1 = _mm256_max_pd(g2, g3);
      __m256d j2 = _mm256_max_pd(g4, g5);
      __m256d j3 = _mm256_max_pd(g6, g7);

      __m256d L0 = _mm256_max_pd(j0, j1);
      __m256d L1 = _mm256_max_pd(j2, j3);

      __m256d L2 = _mm256_max_pd(L0, L1);

      // most compilers optimize this away via register unpacking. 
      // It's the tail end of a horizontal max reduction.
      double q[4];
      _mm256_store_pd(&q[0], L2);
      double q0 = q[0] > q[1] ? q[0] : q[1];
      double q1 = q[2] > q[3] ? q[2] : q[3];
      double q2 = q0 > q1 ? q0 : q1;
      mpr[i] = q2 > mpr[i] ? q2 : mpr[i];
   }
   /*
     These stores serve only 2 purposes.
     First, they are used to compute "edge" components, where
     the range of d is constrained to be < unroll width 
     (we're using 32 for AVX)

     Second, they are used to compute divergence. Normally we explicitly
     compute initial co-moments. One way to test that calculations didn't diverge
     is to explicitly test against the final co-moments.

     This will actually tell you if accumulation failed somewhere.
     Actually it might be worth making that into a standard feature so as to issue an explicit warning.

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

  This is only for constrained cases. 

  Unlike the unconstrained case, we don't rely on
  knowing dcount. 

  rbegin indicates some starting value corresponding
  to an offset.

  rcount corresponds to an upper limit.

  It's assumed that the range of "d" is constrained
  rcount - r, where  0 <= r < rcount

  meaning that only d = 0 executes when r == rcount - 1

  Experience suggests that there isn't an obvious way to make
  optimized implementations appear simple and clear, so
  this focuses on minimizing the number of explicit constraints
  that someone who wishes to use these must interact with.

  Notice that this requires an extra parameter, since we can't assume
  that the maximum range of d is equal to stride - 1 for the stride 
  used in the unconstrained case. For a very short input, it could be less.

*/

static inline void pearson_edge(
    double* __restrict cov,
    double* __restrict mpr,
    double* __restrict mpc,
    double* __restrict dr_bwd,
    double* __restrict dc_bwd,
    double* __restrict dr_fwd,
    double* __restrict dc_fwd,
    double* __restrict invnr,
    double* __restrict invnc,
    int rbegin,
    int rcount){

   int dcount = rcount - rbegin;

   for (int r = 0; r < rcount; ++r){
      for (int d = 0; d < dcount - r; ++d){
         if (r != 0){
            cov[d] -= dr_bwd[r - 1] * dc_bwd[d - 1];
            cov[d] += dr_fwd[r - 1] * dc_fwd[d - 1];
         }
         double cr = cov[d] * invnr[r] * invnc[d];
         if (mpr[r] < cr){
            mpr[r] = cr;
         }
         if (mpc[d] < cr){
            mpc[d] = cr;
         }
      }
   }
}


/*

cv is a sequence of co-moments

mpr is a matrix profile buffer, corresponding to indices 0 ... subsequence count - min separation
mpc is a matrix profile buffer, corresponding to indices min separation .... subsequence count

They should be allocated separately. Aliasing writes can produce weird bugs, and it's undefined behavior when using restrict.

The other difference is that mpr always refers to the subsequence with the lowest index in any comparison. In the initial 
input of cv, it is uniformly 0, and cv naturally has subsequence count - min separation elements as a result, passed here as "dcount",
standing for "diagonal count".


It's assumed that any missing data regions have been trimmed from the beginning and end via preprocessing (see my matlab code).


Formulas for the difference equations:
     "dr_bwd", "dc_bwd", "dr_fwd", "dc_fwd" 

are seen in the Matlab version of this.


They can be rederived pretty easily starting from the covariance update formula at the end of:

Numerically Stable, Scalable Formulas for Parallel and
Online Computation of Higher-Order Multivariate
Central Moments with Arbitrary Weights


*/

void compute_mp(double* __restrict cv,
                double* __restrict mpr,
                double* __restrict mpc,
                double* __restrict dr_bwd,
                double* __restrict dc_bwd,
                double* __restrict dr_fwd,
                double* __restrict dc_fwd,
                double* __restrict invnr,
                double* __restrict invnc,
                int dcount){

   constexpr int stride = 32;

   for(size_t d = 0; d < dcount; d += stride){
      // Compute the number of steps over which we can advance cv[i]...cv[i+stride-1], endpoints included
      int full = dcount - d - stride + 1 > 0 ? dcount - d - stride + 1 : 0;
      if(full != 0){
         pearson_inner(cv + d, mpr, mpc + d, dr_bwd, dc_bwd + d, dr_fwd, dc_fwd + d, invnr, invnc + d, full);
      }
      pearson_edge(cv + d, mpr, mpc + d, dr_bwd, dc_bwd + d, dr_fwd, dc_fwd + d, invnr, invnc + d, full, dcount - d - full);
   }
}


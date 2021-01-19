#include <immintrin.h>
#include <algorithm>
#include <cmath>
#include "mex.h"

/*
 *
 * I am leaving a few notes for colleagues, who may be reading this.
 *
 *
 * First, regarding nans and subsequences containing only constants.
 *
 * The methods at this level have no knowledge of these things.
 * They are optimized for finite data, which has a couple edge cases.
 *
 *
 * Since not everyone is aware of this, I'll point out the source of divergence in cases where intermediate arithmetic produces nan values
 * (typically via overflow followed by multiplication)
 *
 * The most common assembly implementation of max (which is not compliant with IEEE2018 rules)
 *
 * is
 *
 * max(a,b) = a > b ? a : b
 *
 *
 * This means:
 *
 * if a or b are +- nan, then
 *
 * a > b == false
 *
 * and so
 *
 * max(a,b) = b
 *
 * Now if you have a cascaded reduction
 *
 * max(<a0,b0>.....<ak,bk>)
 *
 * it's possible for the largest value to be in a comparison of the form
 *
 * max(largest, nan) --> nan
 *
 * with the next comparison being
 *
 * max(nan, some_other_value)
 *
 * in which case, some_other_value is silently propagated
 *
 * IEEE2018 prefers to propagate nan. Matlab prefers to propagate a finite value, which is error prone.
 *
 *
 * Second note...
 *
 * I have made considerable efforts to improve the stability of lookahead methods wherever possible.
 * I swapped out the method used for difference equations some time ago. This one unfortunately uses 4 arrays rather than 2, but it held up better in stress tests.
 * I provided some preprocessing to trim missing leading and trailign sections and slightly perturb sequences of constants.
 *
 * The last point will be to explicitly detect that divergence has occurred. This is done by explicitly computing the end values for each diagonal
 * excluding the last.
 *
 * This can issue a warning or return the measurements as a "witness" of the calculations.
 *
 * Last note on implementation...
 *
 * This is specialized for AVX2. As it is, this should run quite well on any target that supports AVX2.
 * If I have to include index calculations, I either have to reduce throughput by halving the unroll factor
 * or accept the use of masked writes, which pessimize run times on AMD hardware.
 *
 * I prefer to ignore this as it keeps the implementation simple.
 *
 * Reference
 * https://www.agner.org/optimize/instruction_tables.pdf
 *
 *
 */


// reference version of the main routine, written to match against earlier matlab code for the main routine
//


void mpx_ref(double* mpr, double* mpc, double* cv, double* dr_bwd, double* dc_bwd, double* dr_fwd, double* dc_fwd, double* invnr, double* invnc, int dcount){
    for(int d = 0; d < dcount; ++d){
        double c = cv[d];
        for(int r = 0; r < dcount - d; ++r){
            int k = r + d;
            if(r > 0){
                c -= dr_bwd[r-1] * dc_bwd[k-1];
                c += dr_fwd[r-1] * dc_fwd[k-1];
            }
            double corr = c * invnr[r] * invnc[k];
            if(corr > mpr[r]){
                mpr[r] = corr;
            }
            if(corr > mpc[k]){
                mpc[k] = corr;
            }
        }
    }
}



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
    for (int i = 0; i < iters; ++i)
    {
        if (i != 0)
        {
            __m256d dr_bwd_ = _mm256_broadcast_sd(dr_bwd + i - 1);
            
            // cv - dr_bwd[row-1] * dc_bwd[col-1]
            c0 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i - 1), c0);
            c1 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 3), c1);
            c2 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 7), c2);
            c3 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 11), c3);
            c4 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 15), c4);
            c5 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 19), c5);
            c6 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 23), c6);
            c7 = _mm256_fnmadd_pd(dr_bwd_, _mm256_loadu_pd(dc_bwd + i + 27), c7);
            
            __m256d dr_fwd_ = _mm256_broadcast_sd(dr_fwd + i - 1);
            
            // cv + dr_fwd[row-1] * dc_fwd[col-1]
            c0 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i - 1), c0);
            c1 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 3), c1);
            c2 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 7), c2);
            c3 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 11), c3);
            c4 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 15), c4);
            c5 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 19), c5);
            c6 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 23), c6);
            c7 = _mm256_fmadd_pd(dr_fwd_, _mm256_loadu_pd(dc_fwd + i + 27), c7);
        }
        
        // corr = cv * (1/norm(subseq[row]) * (1/norm(subseq[col]))
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
        
        _mm256_storeu_pd(mpc + i, h0);
        _mm256_storeu_pd(mpc + i + 4, h1);
        _mm256_storeu_pd(mpc + i + 8, h2);
        _mm256_storeu_pd(mpc + i + 12, h3);
        _mm256_storeu_pd(mpc + i + 16, h4);
        _mm256_storeu_pd(mpc + i + 20, h5);
        _mm256_storeu_pd(mpc + i + 24, h6);
        _mm256_storeu_pd(mpc + i + 28, h7);
        
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
     These stores serve 2 purposes.
     First, they are used to compute "edge" components, where
     the range of d is constrained to be < unroll width
     (we're using 32 for AVX)
    
     Second, they are used to compute divergence. Normally we explicitly
     compute initial co-moments. These give the final co-moment values, which are also explicitly computable.
     One way to test whether divergence occurred during calculations
     is to explicitly compute final values and test against these.
    
     This will actually tell you if accumulation failed somewhere.
    
     */
    /*
    _mm256_store_pd(cv, c0);
    _mm256_store_pd(cv + 4, c1);
    _mm256_store_pd(cv + 8, c2);
    _mm256_store_pd(cv + 12, c3);
    _mm256_store_pd(cv + 16, c4);
    _mm256_store_pd(cv + 20, c5);
    _mm256_store_pd(cv + 24, c6);
    _mm256_store_pd(cv + 28, c7);*/
}

/*
 *
 * This next part is only for constrained cases.
 *
 * Unlike the unconstrained case, we don't rely on
 * knowing dcount, (diagonal count, in a hankel function sense).
 *
 * rbegin indicates some starting value corresponding
 * to an offset.
 *
 * rcount corresponds to an upper limit.
 *
 * It's assumed that the range of "d" is constrained
 * rcount - r, where  0 <= r < rcount
 *
 * meaning that only d = 0 executes when r == rcount - 1
 *
 * Experience suggests that there isn't an obvious way to make
 * optimized implementations appear simple and clear, so
 * this focuses on minimizing the number of explicit constraints
 * that someone who wishes to use these must interact with.
 *
 * Notice that this requires an extra parameter, since we can't assume
 * that the maximum range of d is equal to stride - 1 for the stride
 * used in the unconstrained case. For a very short input, it could be less.
 *
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
    
    for (int r = rbegin; r < rcount; ++r){
        for (int d = 0; d < dcount - (r - rbegin); ++d){
            int k = r + d;
            if (r != 0){
                cov[d] -= dr_bwd[r - 1] * dc_bwd[k - 1];
                cov[d] += dr_fwd[r - 1] * dc_fwd[k - 1];
            }
            double cr = cov[d] * invnr[r] * invnc[k];
            if (mpr[r] < cr){
                mpr[r] = cr;
            }
            if (mpc[k] < cr){
                mpc[k] = cr;
            }
        }
    }
}


/*
 *
 * cv is a sequence of co-moments
 *
 * mpr is a matrix profile buffer, corresponding to indices 0 ... subsequence count - min separation
 * mpc is a matrix profile buffer, corresponding to indices min separation .... subsequence count
 *
 * They should be allocated separately. Aliasing writes can produce weird bugs, and it's undefined behavior when using restrict.
 *
 * The other difference is that mpr always refers to the subsequence with the lowest index in any comparison. In the initial
 * input of cv, it is uniformly 0, and cv naturally has subsequence count - min separation elements as a result, passed here as "dcount",
 * standing for "diagonal count".
 *
 *
 * It's assumed that any missing data regions have been trimmed from the beginning and end via preprocessing (see my matlab code).
 *
 *
 * Formulas for the difference equations:
 * "dr_bwd", "dc_bwd", "dr_fwd", "dc_fwd"
 *
 * are seen in the Matlab version of this.
 *
 *
 * They can be rederived pretty easily starting from the covariance update formula at the end of:
 *
 * Numerically Stable, Scalable Formulas for Parallel and
 * Online Computation of Higher-Order Multivariate
 * Central Moments with Arbitrary Weights
 *
 *
 */


void compute_self_mp(double* __restrict cv,
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
    int edgect = dcount % stride;
    int stridect = dcount / stride;
    for(int s = 0; s < stridect; ++s){
        // Compute the number of steps over which we can advance cv[i]...cv[i+stride-1], endpoints included
        int d = s * stride; 
        int full = dcount - d - (stride - 1);
        pearson_inner(cv + d, mpr, mpc + d, dr_bwd, dc_bwd + d, dr_fwd, dc_fwd + d, invnr, invnc + d, full);
        
        pearson_edge(cv + d, mpr, mpc + d, dr_bwd, dc_bwd + d, dr_fwd, dc_fwd + d, invnr, invnc + d, full, dcount - d);
        
    }
    if(edgect != 0){
        int d = stridect * stride;
        pearson_edge(cv + d, mpr, mpc + d, dr_bwd, dc_bwd + d, dr_fwd, dr_fwd + d, invnr, invnc + d, d, dcount - d);
    }
}

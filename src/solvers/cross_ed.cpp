#include<immintrin.h>
#include<cmath>
#include "cross_cov.h"
#ifdef _OPENMP
#include<omp.h>
#endif
const int stride{32};

/*  Normalized euclidean distance computed using the formula
 *  ed[i] = sqrt( ((x[i : i + subseqlen] - mu[i]) / sig[i]  - cmpseq)**2 )
 *
 *  cmpseq must be normalized externally.
 *  Ignoring the cost of the operations themselves, doing it here would cause the compiler to 
 *  either generate more spill code or throttle throughput.
 *
 *  stride/unroll width here was set to be wide enough to use vectorization and hide instruction latency
 *  on a Skylake like architecture, but they should be okay on any x86_64 cpus released 2013 or later made by AMD or Intel.
 *  
 *  This ordering of operations is simple enough to reproduce with decent throughput on ARM using Neon (or Helium with some adjustments), using
 *  the corresponding intrinsics. You can get something close to this using plain C++, with localized arrays or via a polyhedral
 *  compiler using unroll and jam. I found this approach was easier to distribute for now, given that a number of compilers do not yet
 *  fully support OpenMP 4.5 and the ones that do still don't do a good job with sliding window type problems. It seems like they confuse
 *  the SCEV passes or corresponding compilers just don't implement this as an induction pattern, not sure which..
 *
 *  This assumes 2 fmas may be issued per cycle with single use operands issued as memory references.
 *  Numbers may differ slightly on AMD, but this doesn't have anything that is particularly difficult there.
 *
 *  Reference
 *  https://www.agner.org/optimize/instruction_tables.pdf
 *  https://software.intel.com/sites/landingpage/IntrinsicsGuide/
 *
 */

void crossdist(double* __restrict cc, double* __restrict ts, double* __restrict mu, double* __restrict sig, double* __restrict cmpseq, int seqcount, int seqlen){
    int aligned = seqcount - seqcount % stride;
    int count = aligned / stride;
    // This only parallelizes calculations if you use an openmp library.
    // It may slow down a little if the output array isn't aligned to a cache line boundary
    // It would very slow if w
    #pragma omp parallel for
    for(int h = 0; h < count; ++h){
        int i = h * stride; // Older openmp standards require a normalized loop counter, thus this...
        __m256d c0 = _mm256_setzero_pd();
        __m256d c1 = _mm256_setzero_pd();
        __m256d c2 = _mm256_setzero_pd();
        __m256d c3 = _mm256_setzero_pd();
        __m256d c4 = _mm256_setzero_pd();
        __m256d c5 = _mm256_setzero_pd();
        __m256d c6 = _mm256_setzero_pd();
        __m256d c7 = _mm256_setzero_pd();

        for(int j = 0; j < seqlen; ++j){
            int k = i + j;

	    // r[0:stride] = ts[i + j: i + j + stride]
	    __m256d r0 = _mm256_loadu_pd(ts + k);
	    __m256d r1 = _mm256_loadu_pd(ts + k + 4);
            __m256d r2 = _mm256_loadu_pd(ts + k + 8);
            __m256d r3 = _mm256_loadu_pd(ts + k + 12);
            __m256d r4 = _mm256_loadu_pd(ts + k + 16);
            __m256d r5 = _mm256_loadu_pd(ts + k + 20);
            __m256d r6 = _mm256_loadu_pd(ts + k + 24);
            __m256d r7 = _mm256_loadu_pd(ts + k + 28);

            __m256d t0 = _mm256_loadu_pd(mu + i);
            __m256d t1 = _mm256_loadu_pd(mu + i + 4);
            __m256d t2 = _mm256_loadu_pd(mu + i + 8);
            __m256d t3 = _mm256_loadu_pd(mu + i + 12);
            __m256d t4 = _mm256_loadu_pd(mu + i + 16);
            __m256d t5 = _mm256_loadu_pd(mu + i + 20);
            __m256d t6 = _mm256_loadu_pd(mu + i + 24);
            __m256d t7 = _mm256_loadu_pd(mu + i + 28);

            // r[0:stride] = r[0:stride] - mu[i:i+stride]
	    r0 = _mm256_sub_pd(r0, t0);
            r1 = _mm256_sub_pd(r1, t1);
            r2 = _mm256_sub_pd(r2, t2);
            r3 = _mm256_sub_pd(r3, t3);
            r4 = _mm256_sub_pd(r4, t4);
            r5 = _mm256_sub_pd(r5, t5);
            r6 = _mm256_sub_pd(r6, t6);
            r7 = _mm256_sub_pd(r7, t7);

            __m256d q0 = _mm256_loadu_pd(sig + i);
            __m256d q1 = _mm256_loadu_pd(sig + i + 4);
            __m256d q2 = _mm256_loadu_pd(sig + i + 8);
            __m256d q3 = _mm256_loadu_pd(sig + i + 12);
            __m256d q4 = _mm256_loadu_pd(sig + i + 16);
            __m256d q5 = _mm256_loadu_pd(sig + i + 20);
            __m256d q6 = _mm256_loadu_pd(sig + i + 24);
            __m256d q7 = _mm256_loadu_pd(sig + i + 28);

            // r[0:stride] = r[0:stride] / sig[i:i+stride]
	    r0 = _mm256_div_pd(r0, q0);
	    r1 = _mm256_div_pd(r1, q1);
	    r2 = _mm256_div_pd(r2, q2);
	    r3 = _mm256_div_pd(r3, q3);
	    r4 = _mm256_div_pd(r4, q4);
	    r5 = _mm256_div_pd(r5, q5);
	    r6 = _mm256_div_pd(r6, q6);
	    r7 = _mm256_div_pd(r7, q7);

            // r[0:stride] = r[0:stride] - cmpseq[j]
            __m256d s = _mm256_broadcast_sd(cmpseq + j);
 
            r0 = _mm256_sub_pd(r0, s);           
            r1 = _mm256_sub_pd(r1, s);           
            r2 = _mm256_sub_pd(r2, s);           
            r3 = _mm256_sub_pd(r3, s);           
            r4 = _mm256_sub_pd(r4, s);           
            r5 = _mm256_sub_pd(r5, s);           
            r6 = _mm256_sub_pd(r6, s);           
            r7 = _mm256_sub_pd(r7, s);           

	    // ed2[0:stride] = ed2[0:stride] + r[0:stride]**2
	    c0 = _mm256_fmadd_pd(r0, r0, c0);
	    c1 = _mm256_fmadd_pd(r1, r1, c1);
	    c2 = _mm256_fmadd_pd(r2, r2, c2);
	    c3 = _mm256_fmadd_pd(r3, r3, c3);
	    c4 = _mm256_fmadd_pd(r4, r4, c4);
	    c5 = _mm256_fmadd_pd(r5, r5, c5);
	    c6 = _mm256_fmadd_pd(r6, r6, c6);
	    c7 = _mm256_fmadd_pd(r7, r7, c7);
        }

	// Note: square roots are not a typical bottleneck here. Most of the cost of computing this is division
	//       You can optionally pre-compute 1/sig and use multiplication instead. The reference link above should
	//       give you some idea of how ugly division latencies are on modern hardware. Similar to square roots, they
	//       appear to rely on the use of micro coded routines rather than a simple opcode -> micr-op sequence conversion
	

        // ed[0:stride] = sqrt(ed2[0:stride])
        c0 = _mm256_sqrt_pd(c0);
        c1 = _mm256_sqrt_pd(c1);
        c2 = _mm256_sqrt_pd(c2);
        c3 = _mm256_sqrt_pd(c3);
        c4 = _mm256_sqrt_pd(c4);
        c5 = _mm256_sqrt_pd(c5);
        c6 = _mm256_sqrt_pd(c6);
        c7 = _mm256_sqrt_pd(c7);

	// dists[i:i+stride] = ed[0:stride]
	_mm256_storeu_pd(cc + i, c0);
	_mm256_storeu_pd(cc + i + 4, c1);
	_mm256_storeu_pd(cc + i + 8, c2);
	_mm256_storeu_pd(cc + i + 12, c3);
	_mm256_storeu_pd(cc + i + 16, c4);
	_mm256_storeu_pd(cc + i + 20, c5);
	_mm256_storeu_pd(cc + i + 24, c6);
	_mm256_storeu_pd(cc + i + 28, c7);

    }

    // compute the remaining portion, which does align to a multiple of our unroll width
    for(int i = aligned; i < seqcount; ++i){
        double ed = 0;
	for(int j = 0; j < seqlen; ++j){
            int k = i + j;
	    double r = (ts[k] - mu[i])/sig[i] - cmpseq[j];
            ed += r * r;
	}
	cc[i] = sqrt(ed);
    }
}


#include<immintrin.h>
#include "cross_cov.h"

const int stride{32};

void crosscov(double* __restrict cc, double* __restrict ts, double* __restrict mu, double* __restrict cmpseq, int seqcount, int seqlen){
    int aligned = seqcount - seqcount % stride;

    for(int i = 0; i < aligned; i += stride){
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
	    __m256d r0 = _mm256_loadu_pd(ts + k);
	    __m256d r1 = _mm256_loadu_pd(ts + k + 4);
            __m256d r2 = _mm256_loadu_pd(ts + k + 8);
            __m256d r3 = _mm256_loadu_pd(ts + k + 12);
            __m256d r4 = _mm256_loadu_pd(ts + k + 16);
            __m256d r5 = _mm256_loadu_pd(ts + k + 20);
            __m256d r6 = _mm256_loadu_pd(ts + k + 24);
            __m256d r7 = _mm256_loadu_pd(ts + k + 28);

            __m256d t0 = _mm256_loadu_pd(mu + k);
            __m256d t1 = _mm256_loadu_pd(mu + k + 4);
            __m256d t2 = _mm256_loadu_pd(mu + k + 8);
            __m256d t3 = _mm256_loadu_pd(mu + k + 12);
            __m256d t4 = _mm256_loadu_pd(mu + k + 16);
            __m256d t5 = _mm256_loadu_pd(mu + k + 20);
            __m256d t6 = _mm256_loadu_pd(mu + k + 24);
            __m256d t7 = _mm256_loadu_pd(mu + k + 28);

            r0 = _mm256_sub_pd(r0, t0);
            r1 = _mm256_sub_pd(r1, t1);
            r2 = _mm256_sub_pd(r2, t2);
            r3 = _mm256_sub_pd(r3, t3);
            r4 = _mm256_sub_pd(r4, t4);
            r5 = _mm256_sub_pd(r5, t5);
            r6 = _mm256_sub_pd(r6, t6);
            r7 = _mm256_sub_pd(r7, t7);

            __m256d s = _mm256_broadcast_sd(cmpseq + j);
            __m256d c0 = _mm256_fmadd_pd(r0, s, c0);
	    __m256d c1 = _mm256_fmadd_pd(r1, s, c1);
	    __m256d c2 = _mm256_fmadd_pd(r2, s, c2);
	    __m256d c3 = _mm256_fmadd_pd(r3, s, c3);
	    __m256d c4 = _mm256_fmadd_pd(r4, s, c4);
	    __m256d c5 = _mm256_fmadd_pd(r5, s, c5);
	    __m256d c6 = _mm256_fmadd_pd(r6, s, c6);
	    __m256d c7 = _mm256_fmadd_pd(r7, s, c7);
        }

	_mm256_storeu_pd(cc + i, c0);
	_mm256_storeu_pd(cc + i + 4, c1);
	_mm256_storeu_pd(cc + i + 8, c2);
	_mm256_storeu_pd(cc + i + 12, c3);
	_mm256_storeu_pd(cc + i + 16, c4);
	_mm256_storeu_pd(cc + i + 20, c5);
	_mm256_storeu_pd(cc + i + 24, c6);
	_mm256_storeu_pd(cc + i + 28, c7);

    }


    for(int i = aligned; aligned < seqcount; ++aligned){
        double cv = 0;
	for(int j = 0; j < seqlen; ++j){
            int k = i + j;
	    cv += (ts[k] - mu[k]) * cmpseq[j];
	}
	cc[i] = cv;
    }
}




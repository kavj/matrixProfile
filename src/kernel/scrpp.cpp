#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"

// I will clean up the __Restrict__ and other nonsense later. They're mainly used by the compiler during auto-vectorization. 
// There are a lot of different ways to set up this section in order to accommodate extremely long or short subsequence lengths and other things. 

// This interleaves several queries, which turns a sparse problem into a dense one. It allows us to iterate continuously over a memory buffer. 
// It is important to note that we need to be careful about exclusion zones
// Since the queries are in order we naturally hit them in order
// I will probably have to account for that by using a temporary indexed list whenever this is encountered. I would just need to declare a maxunroll or something at compile time and check "if conflicts" 
// This is easy enough to deal with when using simd registers. Just int64_t excluded_mask --> blend if a conflict is detected

// typically unroll is anywhere from 8 to 12 multiplied by intended simd length. It can be modified to ensure it's reasonable with respect to cache. 
// With very long subsequence lengths, it may not be obvious how to tune this. The longest example in the paper was 4096. This should be reasonable up to around 8192 with scalar extensions or 1024 with simd
// It's a reasonable guess at least based on typical L1/L2 sizes used in x86 in recent years
//
//
// It occurs to me that the primary candidates for parallelization would be the overall time series and batch normalization
// We could therefore have a loop for batch norm
// We would of course need some way to partition qcov among threads?
// This probably requires a struct or something somewhere


// This probably ultimately needs some shared data structure like struct ts_descript{ bufferlist .....}
//

// It could also make sense to use an MKL back end? This way we could just set up xcorr directly

void batch_normalize(const double* __Restrict__ ts, const double* __Restrict__ q, const double* __Restrict__ mu, int qstart, int qcount, int sublen, int step){
   int back = len-sublen+1;
   block<double> reg;
   // hoist means 
   for(int i = 0; i < qcount; i++){
      reg(i) = mu[i*step+qstart];
   } 
   for(int i = 0; i < sublen; i++){
      for(int j = 0; j < qcount; j++){
         q[i*qcount+j] = ts[i+qstart+j*step]-reg(j);
      }
   }
}


// it would be easier to track things entering or leaving the exclusion zone since queries are ordered
// Then it becomes prefix ...  partial overlap ... suffix
//
// In general these things are spaced such that we shouldn't have a lot of things in there simultaneously
// In fact I should just fall back to a reference if it's too narrow



static inline void partial_xcorr_kern(double* __Restrict__ qmcov, double* __Restrict__ qcorr, double* __Restrict__ qbuf, int64_t* __Restrict__ qind, const double* __Restrict__ ts, const double* __Restrict__ mu, const double* __Restrict__ invn, int qstart, int unroll, int sublen, int qstart, int step, int excl){
   block<double> cov;
   for(int i = 0; i < len; i++){
      double m = mu[i];
      for(int j = 0; j < sublen; j++){  // beter to avoid checks in tight inner loops, particularly if we use simd later
         double c = ts[i] - m;
         for(int k = 0; k < unroll; k++){
            cov(k) += qbuff[j*unroll+k]*c;
         } 
      }
      // figure out which queries overlap the exclusion zone
      int pref_end = max(0,(std::max(0,i-excl) - qstart)/step);
      int suff_begin = max(0,(std::min(len-sublen+1,i+excl) - qstart)/step);
      for(int k = 0; k < pref_end; k++){
         int qi = qstart+k*step; // <-- I'll have to pass an offset pointer somewhere?
         qcorr[qi] = cov(k);
         qcov[qi] = cov(k-unroll);
         qind[qi] = i;

      }
      for(int k = suff_begin; k < unroll; k++){
         int qi = qstart+k*step; 
         qcorr[qi] = cov(k);
         qcov[qi] = cov(k-unroll);
         qind[qi] = i;

      }
   }
}




// This takes a query range and preallocated buffers
// qcorr could technically be replaced with temporary variables. It's just rescaled from qcov
// In general we only care about the maxima here, which helps cut down on excess bus traffic
void prescr_partial_xcorr(double* __Restrict__ qmcov, double* __Restrict__ qcorr, double* __Restrict__ qbuf, int64_t* __Restrict__ qind, const double* __Restrict__ ts, const double* __Restrict__ mu, const double* __Restrict__ invn, int qstart, int qcount, int len, int sublen, int step, int excl){
   int rem = qcount;
   while(rem > 0){
      int unroll = std::min(8,rem);
      int start = qstart+(qcount-rem)*step;
      batch_normalize(ts,qbuf,mu,start,unroll,sublen,step);
      partial_xcorr_kern(qmcov,qcorr,qbuf,qind,ts,mu,invn,start,step,excl); 
      rem -= unroll;
   }
}




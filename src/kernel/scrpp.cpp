#include "../utils/reg.h"

// I will clean up the __Restrict__ and other nonsense later. They're mainly used by the compiler during auto-vectorization. 
// There are a lot of different ways to set up this section in order to accommodate extremely long or short subsequence lengths and other things. 

// These are very sparse calculations, so simd doesn't help much. It's more down to things like prefetching and grouping. By batching based on unroll factors, we get a chunk of queries into a contiguous chunk of memory. It might be ideal to do this whole function per thread with heavier unrolling


// This would probably be faster if I interleaved queries to allow for simd usage, since we could force aligned loads and broadcast terms of ts
// The copy stage would be a little exotic, but that doesn't really matter. I should probably benchmark compare first?

// Actually these should be grouped into contiguous memory rather than strided. The main problem is that the strided case can produce cache conflicts with long subsequence lengths.
// It's still possible to encounter that problem, but it's much less pronounced if they aren't used again for some time. 
// My suspicion is that I should eventually compare this to a recursive strategy, given that O(m^~1.36) is achievable for m outputs of sublen m.
// This function shouldn't use simd. It should just do a manual gather here.
void batch_normalize(const double* __Restrict__ ts, const double* __Restrict__ q, const double* __Restrict__ mu, const double* __Restrict__ invnm, int qstart, int qcount, int sublen, int step){
   int back = len-sublen+1;
   block<double> reg;
   // hoist the stat variables
   for(int i = 0; i < qcount; i++){
      reg(i) = mu[qstart+i*step];
   } 
   for(int i = 0; i < qcount; i++){
      reg(i+qcount) = invnm[qstart+i*step];
   }
   for(int i = 0; i < sublen; i++){
      for(int j = 0; j < qcount; j++){
         int k = i+j*step;
         q[k] = ts[qstart+k]-reg(j);
      }
      for(int j = 0; j < qcount; j++){
         int k = i+j*step;
         q[k] *= reg(j+qcount);
      }
   }
}


// self similarity or partial auto correlation can be performed using a smaller number of arrays using the skew symmetric formula for df,dg
// For now let's assume centering one is good enough. This will reduce register pressure somewhat since xmm and ymm only have 16 names
//
void prescr_partial_xcorr(double* __Restrict__ qmcov, double* __Restrict__ const double* __Restrict__ queries, ,const double* __Restrict__ ts, const double* __Restrict__ invn, int len, int sublen, int step){
   


}




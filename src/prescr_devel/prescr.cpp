#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"



// This could probably use some cleanup
// centered unit normalization
void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int qstart, int count, int sublen, int step){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qind = i*sublen;
      int tsind = i*step+qstart;
      double qm = mu[qind];
      double qinvn = invn[qind];
      for(int j = 0; j < sublen; j++){
         qbuf[qind] = (ts[tsind]-qm)*qinvn;
         ++qind;
         ++tsind;
      }
   }
}

// Should I actually use this? It may be acceptible to just adjust scale?
void batch_partial_autocov(double* __restrict__ cbuf, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, const int* __restrict__ qind, int qstart, int count, int sublen, int step){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qi = (qstart+i)*step;  
      int ci = qind[qi];
      double cm = mu[ci];
      double qm = mu[qi];
      double cv = 0;
      for(int j = 0; j < sublen; j++){
         cv += (ts[ci] - cm)*(ts[qi] - qm);
         ++qi;
         ++ci;
      }
      cbuf[i] = cv;
   }
}


// prescrimp uses cross correlation, so it's only possible to fuse normalization and comparison
// we have a partially normalized cross correlation. We hoist chunks of 8 loads  then use a tiny decision tree for the intermediate reduction
// I like it because it's geeky
void max_partial_autocorr_reduction(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, int offset, int count, int qbase){
   const int unroll = 8;
   int aligned = count - count%unroll + offset;
   double qcorr = -1.0;
   int qind;
   for(int i = offset; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*invn[i+j];
      }
      for(int j = 0; j < unroll; j++){
         if(xcorr[i+j] < corr(j)){
            xcorr[i+j] = corr(j);
            cindex[i+j] = qbase;
         }
      }
      block<int> cind;
      for(int j = 0; j < 4; j++){
         if(corr(2*j) < corr(2*j+1)){
            corr(2*j) = corr(2*j+1);
            cind(j) = 2*j+1;
         }
         else{
            cind(j) = 2*j;
         }
      }
      for(int j = 0; j < 2; j++){
         if(corr(4*j) < corr(4*j+2)){
            corr(4*j) = corr(4*j+2);
            cind(4*j) = cind(4*j+2);
         }
      }
      if(corr(0) < corr(4)){
         if(qcorr < corr(4)){
            qcorr = corr(4);
            qind = cind(4);
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cind(0);
      }
   }
   for(int i = aligned; i < offset+count; i++){
      double corr = cov[i]*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i;  //<-- this way it's only forward iteration. Doing both consecutively may partially invalidate  hardware prefetching
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbase;
      }
   }
   cindex[qbase] = qind;  // <-- temporary so these aren't optimized out
   xcorr[qbase] = qcorr;
}


// The time complexity of this is linear or close to linear, so we can get it working and ignore
// I changed it to a zigzag pattern. 
//
void maxpearson_extrap_partialauto(const double* __restrict__ qcov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   for(int i = 0; i < count; i++){
      c = qcov[i];
      int qi = qind[i];
      int ci = i*stride;
      int lim = std::min(ci,qi);
      lim = std::min(extraplen,lim);
      double c = qcov[i];
      for(int j = 0; j < lim; j++){
         c -= dx[ci-j]*df[qi-j];
         c -= dx[qi-j]*df[ci-j];
         double scale = invn[qi]*invn[ci];
         double corr = c*scale;
         if(mp[ci] < corr){
            mp[ci] = corr;
            mpi[ci] = qi;
         }
         if(mp[qi] < corr){
            mp[qi] = corr;
            mpi[qi] = ci;
         }
         --ci;
         --qi;
      }
      int lim = std::max(ci,qi);
      lim = std::min(extraplen,len-lim);
      for(int j = 0; j < lim; j++){
         c += dx[ci]*df[qi];
         c += dx[qi]*df[ci];
         double scale = invn[qi]*invn[ci];
         double corr = c*scale;
         if(mp[ci] < corr){
            mp[ci] = corr;
            mpi[ci] = qi;
         }
         if(mp[qi] < corr){
            mp[qi] = corr;
            mpi[qi] = ci;
         }
         ++ci; 
         ++qi;
      }
      
   }
}


// Maybe skip the partitioning. We have strong cache reuse over the query section via df, dg
//


// experimental section
/*typedef __m256d vtype;

void avx_max_corr_reduction(double* __restrict__ cov,  double* __restrict__ xcorr, int* __restrict__ cindex, double* __restrict__ invn, double qcorr, int qind, int qbase, int offset, int count){
   const int simlen = 4;
   const int unroll = 4;
   int aligned = count - count%unroll + offset;
  for(int i = offset; i < aligned; i+= unroll*simlen){
      block<vtype> corr;
     for(int j = 0; j < unroll; j++){
         corr(j) = aload(cov,i+j*simlen)*uload(inv,i+j*simlen);
      }
      block<vtype> mask;
      block<vtype> aux;
      for(int j = 0; j < unroll; j++){
         aux(j) = uload(xcorr,i+j*simlen);
         mask(j) = corr(j) > uload(xcorr,i+j*simlen);
      }
      for(int j = 0; j < unroll; j+=simlen){
         blend(corr(j),aux(j),mask(j));
         aux(j) = blend(
      }
   }
}

// q is typically pre-normalized, but here we need centered only
double autocov_single(const double* __restrict__ ts, double tmu, double qmu, int offset, int qoffset, int sublen){
   double q = 0;
   for(int i = 0; i < sublen; i++){
      q += (ts[offset+i]-tmu)*(ts[qoffset+i]-qmu);
   }
   return q;
}

*/

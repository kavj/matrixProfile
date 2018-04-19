#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"



// This could probably use some cleanup

void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int qstart, int count, int sublen, int step){
   #pragma omp simd  // <-- check whether syntax is valid stacked, possibly inline this without restrict quantifieers
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

void batch_cov(double* __restrict__ cbuf, const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, const int* __restrict__ qind, int qstart, int count, int sublen, int step){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qi = (qstart+i)*step;  
      int ci = qind[qi];
      double cm = mu[ci];
      double qm = mu[qi];
      double cv = 0;
      for(int j = 0; j < sublen; j++){
         cv += (ts[ci] - cm)*(qs[qi] - qm);
         ++qi;
         ++ci;
      }
      cbuf[i] = cv;
   }
}


// prescrimp uses cross correlation, so it's only possible to fuse normalization and comparison
// we have a partially normalized cross correlation. We hoist chunks of 8 loads  then use a tiny decision tree for the intermediate reduction
// I like it because it's geeky
void max_corr_reduction(double* __restrict__ cov,  double* __restrict__ xcorr, int* __restrict__ cindex, double* __restrict__ invn, double qcorr, int qind, int qbase, int offset, int count){
   const int unroll = 8;
   int aligned = count - count%unroll + offset;
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
         qind =  std::max(0,i-sublen);  //<-- this way it's only forward iteration. Doing both consecutively may partially invalidate  hardware prefetching
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbase;
      }
   }
   cindex[qbase] = qind;  // <-- temporary so these aren't optimized out
   xcorr[qbase] = qcorr;
}


// generally easiest to just figure the necessary starting index, then iterate forward only
// put in limits
void maxpearson_extrap(const double* __restrict__ qcov, const double* __restrict__ df, const double* __restrict__ dx, const double* __restrict__ invn, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen){
   for(int i = 0; i < count; i++){
      int cind = qind[i]; 
      int qind = i*stride;
      double c = qcov[i];
      for(int j = 0; j < extraplen; j++){
         c += dx[cind]*df[qind];
         c += dx[qind]*df[cind];
         double scale = (invn[qind]*invn[cind]);
         double corr = c*scale;
         if(mp[cind] < corr){
            mp[qind] = corr;
            mpi[qind] = cind;
         }
         if(mp[qind] < corr){
            mp[cind] = corr;
            mpi[cind] = qind;
         }
         ++cind; 
         ++qindd;
      }
   }
}

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
}*/

// q is typically pre-normalized, but here we need centered only
double autocov_single(const double* __restrict__ ts, double tmu, double qmu, int offset, int qoffset, int sublen){
   double q = 0;
   for(int i = 0; i < sublen; i++){
      q += (ts[offset+i]-tmu)*(ts[qoffset+i]-qmu);
   }
   return q;
}



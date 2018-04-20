#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"



// This could probably use some cleanup
// centered unit normalization
// this is probably an acceptable way to do it, 
// I was debating whether I should loop over singles in a parallel region
// Perhaps should consult elsewhere?


// should be inlined
/*void normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, double mu, double invn, int qstart, int sublen){
   for(int j = 0; j < sublen; j++){
      qbuf[j] = (ts[j]-qm)*qinvn;
   }
}*/


void fast_invcn(double* __restrict__ invn, const double* __restrict__ ts, const double* __restrict__ mu, int len, int sublen){
   double a = 0;
   for(int i = 0; i < sublen; i++){
      double term = ts[i] - mu[0];
      a += term*term;
   }
   invn[0] = 1.0/s;
   for(int i = sublen; i < len-sublen+1; i++){
      double b = ts[i-sublen];
      double c = ts[i];
      s += ((b - mu[i-sublen]) + (c - mu[i]) * (b - c); 
      invn[i-sublen] = 1.0/s;
   }
}


// mmay not be needed
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


// This is probably the best I can do without avx

void max_reduction(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, double qinvn, double qcov, double qcorr, int qind, int qbaseind, int offset, int count){
   const int unroll = 8;
   int aligned = count - count%unroll + offset;
   for(int i = offset; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*invn[i+j];
      }
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*invn[i+j];
      }
      for(int j = 0; j < unroll; j++){
         if(xcorr[i+j] < corr(j)){
            xcorr[i+j] = corr(j);
            cindex[i+j] = qbaseind;
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
            qind = i+cind(4);
            qcov = cov[i+cind(4)];
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cov[i+cind(0)];
      }
   }
   for(int i = aligned; i < offset+count; i++){
      double corr = cov[i]*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i;  //<-- this way it's only forward iteration. Doing both consecutively may partially invalidate  hardware prefetching
         qcov = cov[i];
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbaseind;
      }
   }
   return ; // qcorr, qcov, qind, 
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


/

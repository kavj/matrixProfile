#include<algorithm>
#include<cstdint>
#include<cmath>
#include "../utils/reg.h"
#include "descriptors.h"


// mmay not be needed

void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int count, int sublen, int qstart, int qstride, int qbufstride){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qind = qstart + i*qstride;
      double qm = mu[qind];
      double qinvn = invn[qind];
      for(int j = 0; j < sublen; j++){
         qbuf[i*qbufstride+j] = (ts[qind+j]-qm)*qinvn;
      }
   }
}


/*void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int qstart, int count, int sublen, int step){
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
}*/


// This is probably the best I can do without avx
// I'm not completely happy with this, because it relies on updating a non-const reference

//struct query_stat rescaled_max_reduct(double* __restrict__ cov, double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, int qbaseind, int cindoffset, int count);
struct query_stat rescaled_max_reduct(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, double qinvn, int qbaseind, int count, int cindoffset){
   const int unroll = 8;
   int aligned = count - count%unroll;
   double qcov = NAN;
   double qcorr = -1.0;
   int qind = -1;
   for(int i = 0; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*qinvn;
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
            qind = i+cind(4) + cindoffset;
            qcov = cov[i+cind(4)];
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cov[i+cind(0)] + cindoffset;
         qcov = cov[i+cind(0)];
      }
   }
   for(int i = aligned; i < cindoffset+count; i++){
      double corr = cov[i]*qinvn*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i+cindoffset;  
         qcov = cov[i];
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbaseind;
      }
   }
   struct query_stat s(qcorr,qbaseind);
   return s;
}


// This can probably be merged with the normal version and an unroll size of 1

// need to add structs in here
// This might be reworked later to remove temporary buffers or to wrap the messy formulas in small inline functions, not sure really but it's too hard to read with inline formulas
void maxpearson_ext_auto(const double* __restrict__ qcov, const double* __restrict__ invn, const double* __restrict__ df, const double* __restrict__ dx, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   for(int i = 0; i < count; i++){
      int qi = qind[i];
      int ci = i*stride;
      int lim = std::min(extraplen,std::min(ci,qi));
      double c = qcov[i];
      for(int j = 1; j < lim; j++){
         c -= dx[ci]*df[qi];
         c -= dx[qi]*df[ci];
         double corr = c*invn[qi]*invn[ci];
         if(mp[ci-j] < corr){
            mp[ci-j] = corr;
            mpi[ci-j] = qi-j;
         }
         if(mp[qi] < corr){
            mp[qi-j] = corr;
            mpi[qi-j] = ci-j;
         }
      }
      lim = std::min(extraplen,len-std::max(ci,qi));
      c = qcov[i];
      for(int j = 1; j < lim; j++){
         c += dx[ci]*df[qi];
         c += dx[qi]*df[ci]; 
         double corr = c*invn[qi]*invn[ci];
         if(mp[ci+j] < corr){
            mp[ci+j] = corr;
            mpi[ci+j] = qi+j;
         }
         if(mp[qi+j] < corr){
            mp[qi+j] = corr;
            mpi[qi+j] = ci+j;
         }
      }
   }
}


// if it's not parallel, we allocate fewer buffers assume parallel for interactivity
int maxpearson_partialauto(struct acorr_desc& aux){
   // int iters = aux.required_passes();
   int alias_offset = aux.q.blen-1;  // truncate edge terms of cross correlation
   for(int i = 0; i < iters; i++){     
      int qct = (i == iters - 1) ? aux.q.bcount : aux.tailcount;
      #pragma omp parallel for
      for(int j = 0; j < aux.q.count; j++){
         double* query = aux.q(j);
         int index_st = (i+j)*aux.qbasestride;
         double m = ac.mu[index_st];
         for(int k = 0; k < sublen; k++){
            q0[k] = ac.ts[index_st+k] - m;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j <  ; j++){
         for(int k = 0; k < qct; k++){
            int status = vsldCorrExecX1D(ac.covdesc[j],aux.q(k),1,aux.covbufs(j));
            int qbaseind = (j*aux.q.bcount + k)*aux.qbasestride;
            //int cindoffset =  
            rescaled_max_reduct(aux.qcov(j)+ alias_offset,ac.xcorr+ac.offset(j),ac.invn+offset(j),ac.xind(j),ac.invn[blk(0)+qbaseind],qbaseind,);
         }
      }      
      // reduce over smaller shared buffers here?
   }
   // reduce over thread buffers
   for(int i = 0; i < ; i++){
      
   }
   //maxpearson_extrap_partialauto(qcov,invn,df,dx,qind, mp,mpi,count,stride,extraplen,len);
   return 0;
}



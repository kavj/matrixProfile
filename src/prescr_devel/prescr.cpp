#include<algorithm>
#include<cstdint>
#include "../utils/reg.h"
#include "descriptors.h"


void fast_invcn(double* __restrict__ invn, const double* __restrict__ ts, const double* __restrict__ mu, int len, int sublen){
   double a = 0;
   for(int i = 0; i < sublen; i++){
      double t = ts[i] - mu[0];
      a += t*t;
   }
   invn[0] = 1.0/a;
   for(int i = 1; i < len-sublen+1; i++){
      double b = ts[i+sublen-1];
      double c = ts[i-1];
      a += ((b - mu[i+sublen-1]) + (c - mu[i-1])) * (b - c); 
      invn[i-sublen+1] = sqrt(1.0/a);
   }
}

// mmay not be needed

void batch_normalize(double* __restrict__ qbuf,  const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ invn, int count, int sublen, int qstart, int qstride, int qbufstride){
   #pragma omp parallel for
   for(int i = 0; i < count; i++){
      int qind = qstar + i*qstride;
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
void rescaled_max_reduct(double* __restrict__ cov,  double* __restrict__ xcorr, const double* __restrict__ invn, int* __restrict__ cindex, double &qcov, double& qcorr, int& qind, double qinvn, int qbaseind, int offset, int count){
   const int unroll = 8;
   int aligned = count - count%unroll;
   for(int i = 0; i < aligned; i+= unroll){
      block<double> corr;
      for(int j = 0; j < unroll; j++){
         corr(j) = cov[i+j]*qinvn;;
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
            qind = i+cind(4)+offset;
            qcov = cov[i+cind(4)];
         }
      }
      else if(qcorr < corr(0)){
         qcorr = corr(0);
         qind = cov[i+cind(0)];
      }
   }
   for(int i = aligned; i < offset+count; i++){
      double corr = cov[i]*qinvn*invn[i];
      if(qcorr < corr){
         qcorr = corr;
         qind =  i+offset;  
         qcov = cov[i];
      } 
      if(xcorr[i] < corr){
         xcorr[i] = corr;
         cindex[i] = qbaseind;
      }
   }
}


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
int maxpearson_partialauto(struct p_autocorr& ac, struct corr_auxbuf& aux, int dlen, int qstride){
   int rem = aux.qcount;
   int iters = a aux.querycount/aux.q.count; 
   for(int i = 0; i < iters; i++){     
      int qct = (i == iters - 1) ? aux.q.count : aux.querycount  - (iters-1)*aux.q.count;
      #pragma omp parallel for
      for(int j = 0; j < aux.q.count; j++){
         int k = (i+j)*querystride;
         int m = mu[k];
         double* p = aux.q(j);
         for(int r = 0; r < sublen; r++){
            p[r] = ts[k+r] - m;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j < qs.blkct; j++){
         for(int k = 0; k < qct; k++){
            int status = vsldCorrExecX1D(acd.covdesc[j],qb.q[k*qb.blkstrd],1,);
            // need debugging info on failure
            rescaled_max_reduct(pcorrbuf.cov+j,acd.xcorr,invn,cindex,qcov, qcorr, qind, qinvn, qbaseind, offset = k, count);
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




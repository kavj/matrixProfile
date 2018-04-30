#include<algorithm>
#include<cstdint>
#include<cmath>
#include "../utils/reg.h"
#include "descriptors.h"
#include "../kernel/prescr_Pearson.h"



// if it's not parallel, we allocate fewer buffers assume parallel for interactivity
int maxpearson_partialauto(struct acorr_desc<double>& aux){
   int iters = aux.required_passes();
   //int alias_offset = aux.sublen-1;  // truncate edge terms of cross correlation
   int sublen = aux.sublen;
   for(int i = 0; i < iters; i++){     
      int qct = (i == iters - 1) ? aux.q.bcount : aux.tailqcount;
      #pragma omp parallel for
      for(int j = 0; j < qct; j++){
         double* query = aux.q(j);
         double m = (aux.mu(0)[(i+j)*aux.qbasestride]);
         //double* a = aux.ts(i+j,aux.qbasestride);
         for(int k = 0; k < sublen; k++){
        //    query[k] = a[k] - m;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j <  aux.bcount; j++){
         int cct = (j < aux.bcount-1) ? aux.blen : aux.taillen;
         int cindoffset = j*aux.bstride; 
         for(int k = 0; k < qct; k++){
            int status = vsldCorrExecX1D(aux.covtsks[j],aux.q(k),1,aux.covbufs(j),aux.blen+sublen-1);
            int qbaseind = (j*aux.q.bcount + k)*aux.qbasestride;
            // We should ensure that in the case where this is built as a single threaded operation, we don't allocate multiple buffers and simply reduce queries in place
            rescaled_max_reduct(aux.cov(j),aux.xcorr(j),aux.invn(j),aux.xind(j),aux.invn(0)[qbaseind],qbaseind,aux.blen-sublen+1,j*aux.blen); //need count, need cindoffset);
         }
      }      
      // reduce over smaller shared buffers here?
      for(int j = 0; j < qct; j++){
         // reduce over number of per thread/section buffers
      }
   }
   // reduce over thread buffers
   for(int i = 0; i < 0; i++){
    // maxpearson_ext_auto(aux.q(      
//void maxpearson_ext_auto(const double* __restrict__ qcov, const double* __restrict__ invn, const double* __restrict__ df, const double* __restrict__ dx, const int* __restrict__ qind, double* __restrict__ mp, int* __restrict__ mpi, int count, int stride, int extraplen, int len){
   }
   return 0;
}



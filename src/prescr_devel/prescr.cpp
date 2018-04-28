#include<algorithm>
#include<cstdint>
#include<cmath>
#include "../utils/reg.h"
#include "descriptors.h"
#include "../kernel/prescr_Pearson.h"



// if it's not parallel, we allocate fewer buffers assume parallel for interactivity
int maxpearson_partialauto(struct acorr_desc<double>& aux){
    int iters = aux.required_passes();
   int alias_offset = aux.sublen-1;  // truncate edge terms of cross correlation
   for(int i = 0; i < iters; i++){     
      int qct = (i == iters - 1) ? aux.q.bcount : aux.tailqcount;
      #pragma omp parallel for
      for(int j = 0; j < aux.q.count; j++){
         double* query = aux.q(j);
         double m = *(ac.smu(i+j,ac.qbasestride));
         double* a = ac.ts(i+j,ac.qbasestride);
         for(int k = 0; k < sublen; k++){
            query[k] = a[k] - m;
         }
      }
      #pragma omp parallel for
      for(int j = 0; j <  ac.bcount; j++){
         int cct = (j < ac.bcount-1) ? ac.blen : taillen;
         int cindoffset = j*bstride; 
         for(int k = 0; k < qct; k++){
            int status = vsldCorrExecX1D(ac.covdesc[j],aux.q(k),1,aux.covbufs(j));
            int qbaseind = (j*aux.q.bcount + k)*aux.qbasestride;
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



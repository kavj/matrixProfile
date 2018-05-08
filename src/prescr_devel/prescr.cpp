#include<algorithm>
#include<cstdint>
#include<cmath>
#include "../utils/reg.h"
#include "descriptors.h"
#include "prescr_Pearson.h"


// in general we can bound the amount of memory needed for buffering queries if we aren't using some ridiculously short stride
// in that case it makes more sense to use exact calculation anyway
//

#define buflen 65536


template<typename dtype, typename itype>
int prescrpearson_partialauto(const stridedbuf<dtype>& ts, int sublen){

   const int qstride = sublen/4;
   const int qct = (ts.len-sublen+1)/qstride;  // using this a default base stride
   const int qbuflen = paddedlen(sublen,32);    
   //const int buflen = 
   const int bufct = ts.len/buflen + (ts.len%buflen ? 1 : 0);
   const int basestride = buflen - sublen + 1;
   stridedbuf<dtype> qbuf(qbuflen,qct);
   stridedbuf<dtype> mu(mlen,ts.stride,);
   stridedbuf<dtype> invn(mlen,ts.stride,);
   stridedbuf<dtype> covbufs(buflen,);
   stridedbuf<dtype> qcov(  )
   stridedbuf<dtype> qcorr(  )
   stridedbuf<dtype> qind(  )
   
   //stridedbuf<dtype> mp(mlen,ts.stride,);
   //stridedbuf<itype> mpi(mlen,ts.stride,);
 
   //need the MKL task pointers
  // stridedbuf<>;
   
   xmean_windowed(ts(0),mu(0),ts.len,sublen); 
  // init_dfdx(ts(0),mu(0),df(0),dx(0),sublen,ts.len);
   
   #pragma omp parallel for
   for(int j = 0; j < qnum; j++){
      int skipcount = qbufct*i+j;
      center_query(ts(skipcount,qstride),qbuf(j),mu(skipcount),sublen,qbufct,qstride,qbuflen);
   }
   #pragma omp parallel for 
   for(int i = 0; i < ts.bcount; i++){           
      for(int j = 0; j <  qct; j++){
         int status = vsldCorrExecX1D(covtsks[j],qbuf(j),1,covbufs(j),ts.stride+sublen-1);
         for(int k = 0; k < qbufct; k++){
            
         }
      }
      //int cindoffset = j*aux.bstride; 
      for(int k = 0; k < qct; k++){
          // int status = vsldCorrExecX1D(aux.covtsks[j],aux.q(k),1,aux.covbufs(j),aux.blen+sublen-1);
         //int qbaseind = (j*aux.q.bcount + k)*aux.qbasestride;
         //rescaled_max_reduct(aux.covbufs(j),aux.xcorr(j),aux.invn(j),aux.xind(j),aux.invn(0)[qbaseind],qbaseind,aux.blen-sublen+1,j*aux.blen); //need count, need cindoffset);
      }
   }
   // reduce over thread buffers      
   for(int i = 0; i < qct; i++){

   }
  // stridedbuf<dtype> df(mlen,ts.stride,);
  // stridedbuf<dtype> dx(mlen,ts.stride,);
   
   for(int i = 0; i < qct; i++){
      // extrapolate pass
   }

   return 0;
}


template<typename dtype, typename itype>
int prescrpearson_partialcross(){
   
}
// Add pearson cross correlation here

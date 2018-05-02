#pragma once
#define _POSIX_C_SOURCE 200809L
#include<cstdlib>
#include<mkl.h>
// update this later to deal with further edge cases
// this is only used for the aligned allocations. Use regular malloc for outputs


void* init_buffer(int buflen, int bufct, int alignmt){
   void* buf;
   posix_memalign(&buf,buflen*bufct,alignmt);
   return (void*)buf;
}

// I need to handle error propagation... For now though, it would be good to 
// this should return more than a pointer considering that the calculations basically belong here?
// otherwise we can just pass in everything required I guess

VSLCorrTaskPtr* init_taskptrs(const double* ts, int baselen, int sublen, int blkct, int taillen){
  
    // this looks weird, see MKL docs later 
   VSLCorrTaskPtr* covtsks = (VSLCorrTaskPtr*) malloc(sizeof(VSLCorrTaskPtr)*blkct);
 //  VSLCorrTaskPtr* covtsks = init_buffer(sizeof(VSLCorrTaskPtr),blkct,sizeof(VSLCorrTaskPtr));
   if(covtsks == NULL){
      return nullptr;
   }
   for(int i = 0; i < blkct; i++){
      int status;
      if(i < blkct - 1){
         status = vsldCorrNewTaskX1D(covtsks+i,VSL_CORR_MODE_FFT,baselen,sublen,baselen+sublen-1,ts+i*(baselen-sublen+1),1);
      }
      else{
         status = vsldCorrNewTaskX1D(covtsks+i,VSL_CORR_MODE_FFT,taillen,sublen,taillen+sublen-1,ts+i*(baselen-sublen+1),1);
      }
      if(status != VSL_STATUS_OK){
         // Build up real error checking once we have a less ad hoc architecture in place
//         perror("could not initialize tasks");
         return nullptr; 
      }
   }
   return covtsks;
}

int dest_taskptrs(VSLCorrTaskPtr* covtsks, int blkct){
   for(int i = 0; i < blkct; i++){
     //check this part
     // int status = vsldCorr
     //
   }
   return 0;
}





// update this later to deal with further edge cases
// this is only used for the aligned allocations. Use regular malloc for outputs


void* init_buffer(int bufsz, int alignmt){
   void* buf;
   posix_memalign(buf,buflen*bufct,alignmt);
   return buf;
}

// I'll eventually build in a way to propagate errors from low level libraries to high level bindings. For now I'm mostly ignoring it.
//
VSLCorrTaskPtr* init_taskptrs(const double* ts, int len, int blkct,int blkstride,int sublen){
   VSLCorrTaskPtr* covtsks = init_buffer(1,blockct*sizeof(VSLCorrTaskPtr)*blkct,sizeof(VSLCorrTaskPtr));
   for(int i = 0; i < blockct; i++){
      int status = vsldCorrNewTaskX1D(covtsks+i,VSL_CORR_MODE_FFT,blocklen,sublen,blocklen,ts+i*stride,1);
      if(status != VSL_STATUS_OK){
         // Build up real error checking once we have a less ad hoc architecture in place
         perror("could not initialize tasks");
         return -1; 
      }
   }
   return covtsks;
}

int dest_taskptrs(VSLCorrTaskPtr* covtsks, int blkct){
   for(int i = 0; i < blockct; i++){
     //check this part
     // int status = vsldCorr
     //
   }
   return 0;
}



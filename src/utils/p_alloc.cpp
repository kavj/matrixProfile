#define _POSIX_C_SOURCE 200809L
#include<cstdlib>
#include<iostream>

void* alloc_aligned_buffer(long long buflen, long long alignmt){
   void* buf;
   int chk = posix_memalign((void**)&buf, alignmt, buflen);
   if(chk != 0){  // this could still include results from errno?
      std::cerr << "error allocating memory" << std::endl;
      exit(1);
   }
   return buf;
}

void dealloc_aligned_buffer(void* buf){
   free(buf);
}

// need to find a place for this silly glorified utility macro. 
int paddedlen(long long buflen, long long alignmt){
   return buflen + (buflen % alignmt ? alignmt - buflen % alignmt : 0);
}


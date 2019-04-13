#include "defs.h"
#include "alloc.h"

size_t paddedLen(size_t alignment, size_t byteLen){
   return byteLen + (byteLen % alignment ? alignment - byteLen % alignment : 0);
}

#if CAN_USE_ALIGNED_MALLOC
// windows
#include<malloc.h>

#else

#include<cstdlib>

#endif


void* allocMem(size_t alignment, size_t byteLen){
   #if CAN_USE_ALIGNED_MALLOC
  
   return _aligned_malloc(byteLen, alignment);

   #elif CAN_USE_ALIGNED_ALLOC

   return aligned_alloc(alignment, byteLen);

   #elif CAN_USE_POSIX_MEMALIGN 

   void* dat;
   int err = posix_memalign((void**)&dat, alignment, byteLen);
   return err == 0 ? dat : nullptr;; 

   #else
   
   return malloc(byteLen);

   #endif
}


void deallocMem(void* dat){
#if CAN_USE_ALIGNED_MALLOC
   _aligned_free(dat);
#else
   free(dat);
#endif
}


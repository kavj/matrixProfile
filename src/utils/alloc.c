#include<stdio.h>
#include<errno.h>
#include "alloc.h"

int paddedlen(int buflen){
   return buflen + (buflen % prefalignmt ? prefalignmt - buflen % prefalignmt : 0);
}
// Todo: These could just use malloc wherever 16 byte alignment is sufficient. 

#ifdef _WIN32
#include<stdlib.h>

void* alloc_buff(int buflen){
   void* buf = _aligned_malloc(buflen, prefalignmt);
   if(buf == NULL){
      fprintf(stderr,"problem allocating memory\n");	   
      exit(1);
   }
   return buf;
}

void dealloc_buff(void* v){
   _aligned_free(v);
}

#else

#define _POSIX_C_SOURCE 200809L
#include<stdlib.h>

void* alloc_buff(int buflen){
   void* buf;
   int chk = posix_memalign((void**)&buf, prefalignmt, buflen);
   if(chk != 0){
      if(buf == NULL){
         if((chk == EINVAL) || (chk == ENOMEM)){
            perror("posix_memalign");
            exit(1);
	 }
      }
   }
   return buf;
}

void dealloc_buff(void* v){
   free(v);
}

#endif

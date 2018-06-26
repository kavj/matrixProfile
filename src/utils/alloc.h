


#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__
#define _POSIX_C_SOURCE 200809L
#include<cstdlib>
#include<cstdio>
#include<cerrno>

//Todo: This should be generalized 

inline void* init_buffer(int buflen, int alignmt){
   void* buf;
   int chk = posix_memalign((void**)&buf,alignmt,buflen);
   if(chk != 0){
      if((chk == EINVAL) || (chk == ENOMEM)){
         perror("posix_memalign");
         exit(1);
      }
   }
   return buf;
}


#endif

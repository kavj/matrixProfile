


#define _POSIX_C_SOURCE 200809L
#include<cstdlib>
#include<cstdio>
#include<cerrno>

//Todo: This should be generalized 
// Also: Marked these inline to avoid multiple inclusions due to inlining of other functions. Find a better solution for this.
 
int paddedlen(int buflen, int alignmt){
   return buflen + (buflen%alignmt ? alignmt - buflen%alignmt : 0);
}

void* init_buffer(int buflen, int alignmt){
   void* buf;
   int chk = posix_memalign((void**)&buf,alignmt,buflen);
   if(chk != 0){
      if((chk == EINVAL) || (chk == ENOMEM)){
         perror("posix_memalign");
         exit(1);
      }
      else{
         printf("unknown error\n");
         exit(1);
      }
   }
   return buf;
}

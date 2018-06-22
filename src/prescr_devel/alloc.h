


#ifndef __ALLOC_FILE__
#define __ALLOC_FILE__
#define _POSIX_C_SOURCE 200809L
#include<cstdlib>
#include<cstdio>



inline void* init_buffer(int buflen, int alignmt){
   void* buf;
   int chk = posix_memalign((void**)(&buf),alignmt,buflen);
   //int chk = posix_memalign((void**)(&buf),64,131072*sizeof(double));
   //int chk = posix_memalign((void**)(&buf),alignmt,buflen);
   if(chk != 0){
      //if(chk == EINVAL){
        // printf("unaligned\n");
      //}
      //else if(chck == ENOMEM){
      
       perror("posix_memalign");

     }
  /// }
   //double* b =;
//   exit(1);
   return buf;
}


#endif

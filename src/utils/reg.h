
#ifndef __BLOCK__
#define __BLOCK__
#include <immintrin.h>
typedef __m256d dt;

struct fblock{
   __m256d  r0;
   __m256d  r1;
   __m256d  r2;
   __m256d  r3;
   __m256d  r4;
   __m256d  r5;
   __m256d  r6;
   __m256d  r7;
   __m256d  r8;
   __m256d  r9;
   __m256d  r10;
   __m256d  r11;
   __m256d  r12;
   __m256d  r13;
   __m256d  r14;
   __m256d  r15;
 
   inline __m256d& operator()(int i) __attribute__((always_inline)){
      switch(i){
        case 0:
           return r0;
        case 1:
           return r1;
        case 2: 
           return r2;
        case 3:
           return r3;
        case 4:
           return r4;
        case 5:
           return r5;
        case 6:
           return r6;
        case 7:
           return r7;
        case 8:
           return r8;
        case 9:
           return r9;
        case 10:
           return r10;
        case 11:
           return r11;
        case 12:
           return r12;
        case 13:
           return r13;
        case 14:
           return r14;
        case 15:
           return r15;
      }
   } 
};


struct iblock{
   __m256i r0;
   __m256i r1;
   __m256i r2;
   __m256i r3;
   __m256i r4;
   __m256i r5;
   __m256i r6;
   __m256i r7;
   __m256i r8;
   __m256i r9;
   __m256i r10;
   __m256i r11;
   __m256i r12;
   __m256i r13;
   __m256i r14;
   __m256i r15;
 
   inline __m256i& operator()(int i) __attribute__((always_inline)){
      switch(i){
        case 0:
           return r0;
        case 1:
           return r1;
        case 2: 
           return r2;
        case 3:
           return r3;
        case 4:
           return r4;
        case 5:
           return r5;
        case 6:
           return r6;
        case 7:
           return r7;
        case 8:
           return r8;
        case 9:
           return r9;
        case 10:
           return r10;
        case 11:
           return r11;
        case 12:
           return r12;
        case 13:
           return r13;
        case 14:
           return r14;
        case 15:
           return r15;
      }
   } 
};

#endif

#include<immintrin.h>
#include "../arch/v256.h"


template<typename dt>
struct reg_block{
  // inline reg_block() __attribute__((always_inline)){
  //  }
   dt r0;
   dt r1;
   dt r2;
   dt r3;
   dt r4;
   dt r5;
   dt r6;
   dt r7;
   dt r8;
   dt r9;
   dt r10;
   dt r11;
   dt r12;
   dt r13;
   dt r14;
   dt r15;
 
   inline dt& operator()(int i) __attribute__((always_inline)){
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





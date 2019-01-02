#include "alloc.h"
#include<cstdlib>
#include<iostream>

const int prefalign = 64;

int paddedlen(int buflen){
   return buflen + (buflen % prefalign ? prefalign - buflen % prefalign : 0);
}

void* alloc_aligned_buffer(int buflen){
   void* buf = aligned_alloc(prefalign, buflen);
   if(buf == nullptr){
      std::cerr << "unable to allocate memory" << std::endl;
      exit(1);
   }
   return buf;
}

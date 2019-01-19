#ifndef STRIDEDBUF
#define STRIDEDBUF
#include <algorithm>
#include "alloc.h"


template<typename dtype> struct tilebuf{
   dtype* dat;
   int len;

   tilebuf(int L) : len(L){
      int padlen = paddedlen(L, prefalign);
      dat = reinterpret_cast<dtype*>(alloc_aligned_buffer(padlen * sizeof(dtype), prefalign));
   }

   tilebuf(int L, dtype fillval) : len(L) {
      int padlen = paddedlen(len, prefalign);
      dat = reinterpret_cast<dtype*>(alloc_aligned_buffer(padlen * sizeof(dtype), prefalign));
      if(dat != nullptr){
         std::fill(dat, dat + len, fillval); 
      }
   }

   inline __attribute__((always_inline)) ~tilebuf(){
      dealloc_aligned_buffer(static_cast<void*>(dat));
   }

   inline __attribute__((always_inline)) bool valid() const {
      return (dat != nullptr);
   }

   inline dtype* __attribute__((always_inline)) operator()(){
      return dat;
   }
   
   inline dtype* __attribute__((always_inline)) operator()(int offst){
      return offst < len ? dat + offst : nullptr;
   }

};

typedef tilebuf<double> bufd;
typedef tilebuf<int> bufi;

#endif

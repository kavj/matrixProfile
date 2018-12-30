#ifndef STRIDEDBUF
#define STRIDEDBUF
#include <algorithm>
#include "alloc.h"


template<typename dtype> struct multibuf{
   dtype* dat;
   int bcount;
   int stride;

   multibuf(int count, int blocklen) : bcount(count) {
      stride = paddedlen(blocklen, prefalign);
      dat = reinterpret_cast<dtype*>(alloc_aligned_buffer(stride * count * sizeof(dtype), prefalign));
   }
   
   inline __attribute__((always_inline)) ~multibuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) bool valid() const {
      return (dat != nullptr);
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return i < bcount ? dat + i * stride : nullptr;
   }

};


template<typename dtype> struct primbuf{
   //Todo: allow for externally allocated memory and OS specific allocators

   dtype* dat;
   int len;

   primbuf(int buflen, int fillval) : len(buflen) {
      dat = (buflen > 0) ? reinterpret_cast<dtype*>(alloc_aligned_buffer(paddedlen(buflen, prefalign) * sizeof(dtype), prefalign)) : nullptr;
      if(dat != nullptr){
         std::fill(dat, dat + len, fillval); 
      }
   }

   primbuf(int buflen) : len(buflen) {
      dat = (buflen > 0) ? reinterpret_cast<dtype*>(alloc_aligned_buffer(paddedlen(buflen, prefalign) * sizeof(dtype), prefalign)) : nullptr;
   }
   
   inline __attribute__((always_inline)) ~primbuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) bool valid() const {
      return (dat != nullptr) && (len > 0);
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return (i < len) ? (dat + i) : nullptr;
   }

   // overload of long long is necessary for cases where we need to match the size in bytes of an indexing type and a floating point data type 
   inline dtype*  __attribute__((always_inline)) operator()(long long i){ 
     return (i < len) ? (dat + i) : nullptr;
   }
 
};

#endif

#ifndef STRIDEDBUF
#define STRIDEDBUF
#include <algorithm>
#include "alloc.h"


// We really only need one of these. Just make it default to a unit stride. Stride and len should probably be private with accessors.
// Even if this needed to stop being inline everywhere, it would be fine as this is used to pass pointers to lower level functions anyway

// Note to self, this should be updated for compatibility with msvc
template<typename dtype> struct strided_buffer{
   dtype* dat;
   const int len;
   int stride;

   strided_buffer(int len, int stride) : len(len), stride(stride) {
      dat = reinterpret_cast<dtype*>(alloc_aligned_buffer(len * sizeof(dtype)));
   }

   strided_buffer(int len, int stride, dtype fillval) : len(len), stride(stride) {
      dat = reinterpret_cast<dtype*>(alloc_aligned_buffer(len * sizeof(dtype)));
      if(dat != nullptr){
         std::fill(dat, dat + len, fillval); 
      }
   }

   inline __attribute__((always_inline)) ~strided_buffer(){
      dealloc_aligned_buffer(dat);
   }
  
   inline void __attribute__((always_inline)) set_stride(int s){
      stride = s;
   }

   inline __attribute__((always_inline)) bool valid() const {
      return (dat != nullptr);
   }
 
   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return i * stride < len ? dat + i * stride : nullptr;
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i, int j){ 
      return (i * stride  + j < len) ? dat + i * stride + j : nullptr;
   }

   inline dtype*  __attribute__((always_inline)) operator()(long long i){ 
     return i * stride < len ? (dat + i) : nullptr;
   }
   // test whether this could just be replaced with size_t without incurring a bunch of convert to nonsense, as that tends to break loop vectorization
   inline dtype*  __attribute__((always_inline)) operator()(long long i, long long j){ 
     return (i * stride + j < len) ? dat + i * stride + j : nullptr;
   }

};

#endif

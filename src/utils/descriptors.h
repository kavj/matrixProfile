
#ifndef STRIDEDBUF
#define STRIDEDBUF
#include "alloc.h"
#define prefalign 64


// Note to self: can't use exceptions on object creation, because even if they are correctly propagated, this might result in memory leaks

// rearrange attributes
static inline int   __attribute__((always_inline)) paddedlen(int buflen, int alignmt){return buflen + (buflen%alignmt ? alignmt - buflen%alignmt : 0);}

//I would prefer to avoid the use of late binding object types here due to the use of templating to overload precision of datatypes. This will probably be sufficient.


template<typename dtype> struct multibuf{

   dtype* dat;
   int bcount;
   int stride;

   multibuf(int count, int blocklen) : bcount(count) {
      stride = paddedlen(blocklen, prefalign);
      dat = reinterpret_cast<dtype*>(init_buffer(stride * count * sizeof(dtype), prefalign));
   }
   
   inline __attribute__((always_inline)) ~multibuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) bool isvalid() const {
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

   primbuf(int buflen) : len(buflen) {
      dat = reinterpret_cast<dtype*>(init_buffer(paddedlen(buflen, prefalign) * sizeof(dtype), prefalign));
   }
   
   inline __attribute__((always_inline)) ~primbuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) bool isvalid() const {
      return (dat != nullptr) && (len > 0);
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return (i < len) ? (dat + i) : nullptr;
   }
   
};

#endif

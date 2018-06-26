
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
      stride = paddedlen(blocklen,prefalign);
      dat = reinterpret_cast<dtype*>(init_buffer(stride*count*sizeof(dtype),prefalign));
   }
   
   inline __attribute__((always_inline)) ~multibuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) bool isvalid() const {
      return (dat != nullptr);
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return i < bcount ? dat + i*stride : nullptr;
   }

};


template<typename dtype> struct stridedbuf{
   // This does allow for partial aliasing, which should obviously only be used in the case of pasing read only data structures to different threads
   // If it's externally allocated, then we just want to let it determine the max number of blocks, we can do some minor error checking anyway
   dtype* dat;
   int bcount;
   int stride;
   int len;

   // stridedbuf(dtype* a, int len, int stride, int count) : dat(a), len(len), stride(stride), bcount(count) {} 

   stridedbuf(int buflen) : stride(buflen), len(buflen), bcount(1) {
   //   dat = (dtype*)init_buffer(paddedlen(buflen,prefalign)*sizeof(dtype),prefalign);
      dat = reinterpret_cast<dtype*>(init_buffer(paddedlen(buflen,prefalign)*sizeof(dtype),prefalign));
   }
   
   inline __attribute__((always_inline)) ~stridedbuf(){
      free(dat);
   }

   inline __attribute__((always_inline)) void setstride(int st){
      stride = st;
      bcount = len/st + (len%stride ? 1 : 0);
   }

   inline __attribute__((always_inline)) bool isvalid() const {
      return (dat != nullptr) && (stride*bcount <= len);
   }

   inline __attribute__((always_inline)) int taillen() const {
      return len%stride;
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return (i < bcount) ? (dat + i*stride) : nullptr;
   }
   
};

#endif

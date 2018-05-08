
#ifndef STRIDEDBUF
#define STRIDEDBUF
#include "alloc.h"
#define prefalign 64


// Note to self: can't use exceptions on object creation, because even if they are correctly propagated, this might result in memory leaks

// rearrange attributes
static inline int   __attribute__((always_inline)) paddedlen(int buflen, int alignmt){return buflen + (buflen%alignmt ? alignmt - buflen%alignmt : 0);}


//I would prefer to avoid the use of late binding object types here due to the use of templating to overload precision of datatypes. This will probably be sufficient.

template<typename dtype> struct stridedbuf{
   //stridedbuf(int blen, int bcount, int stride){}
   dtype* dat;
   int bcount;
   int len;
   int stride;
   // This does allow for partial aliasing, which should obviously only be used in the case of pasing read only data structures to different threads
   // If it's externally allocated, then we just want to let it determine the max number of blocks, we can do some minor error checking anyway
   stridedbuf(dtype* a, int len, int stride, int count) : dat(a), len(len), stride(stride), bcount(count), ownsmem(false) {} 

   stridedbuf(int count, int stride) : len(count*stride), bcount(count), stride(stride), ownsmem(true), dat(reinterpret_cast<dtype*>(init_buffer(paddedlen(sizeof(dtype)*count*stride,prefalign),prefalign))) {} 

   stridedbuf(int count, int stride, int len) : len(count*stride), bcount(count), stride(stride), ownsmem(true), dat(reinterpret_cast<dtype*>(init_buffer(paddedlen(sizeof(dtype)*count*stride,prefalign),prefalign))) {}

   inline __attribute__((always_inline)) ~stridedbuf(){
      if(ownsmem){  //set this to something a bit more general 
         free(dat);
      }
   }

   inline __attribute__((always_inline)) bool isvalid() const {
      return (dat != nullptr) && (stride*bcount <= len);
   }

   inline __attribute__((always_inline)) int taillen() const {
      return len%stride;
   }

   inline dtype*  __attribute__((always_inline)) operator()(int i){ 
      return i < bcount ? dat + i*stride : nullptr;
   }

   inline dtype* __attribute__((always_inline)) operator()(int i, int stride){
      return i*stride < len ? dat + i*stride : nullptr;
   }

   private:
   bool ownsmem;
};



#endif


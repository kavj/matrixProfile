
#ifndef PACKED_BUFFER
#define PACKED_BUFFER
#include<new>
#include "alloc.h"


template<typename T> 
class basicBuffer final{
public:
   basicBuffer(std::ptrdiff_t len) : dat_(alloc(len)), len_(len), ownsmem_(true){}
   basicBuffer(T* dat, std::ptrdiff_t len) : dat_(dat), len_(len), ownsmem_(false){} 
   basicBuffer(basicBuffer& bb) = delete;
   basicBuffer operator = (basicBuffer& bb) = delete;
   basicBuffer(basicBuffer&& bb) = default;
      
   ~basicBuffer(){
      if(ownsmem_){
         deallocMem(dat_);
      }
   }
   using iterator = T*;
   using const_iterator = const T*;
   iterator begin(){
      return iterator(dat_);
   }
   iterator end(){
      return iterator(dat_ + len_);
   }
   std::ptrdiff_t len() const{ 
      return len_; 
   }
   bool valid() const{
      return len_ > 0 && dat_ != nullptr;
   }
   T* yieldMem(){  // yield ownership of the underlying buffer
      if(ownsmem_){
         ownsmem_ = false;
         return dat_;
      }
      return nullptr;
   }
   basicBuffer deepCopy(){
      basicBuffer bb(len_);
      if(bb.valid()){
         std::copy(bb.dat_, bb.dat_ + len_, dat_);
      }
      return bb;
   }
   basicBuffer lower(std::ptrdiff_t elemtCount, bool deepCopy = false){
      if(elemtCount > len_){
         return basicBuffer(nullptr, 0);
      }
      basicBuffer bb(dat_, elemtCount);
      return deepCopy ? bb.deepCopy() : bb;
   }
   basicBuffer upper(std::ptrdiff_t elemtCount, bool deepCopy = false){
      if(elemtCount > len_){
         return basicBuffer(nullptr, 0);
      }
      basicBuffer bb(dat_ + (len_ - elemtCount), elemtCount);
      return deepCopy ? bb.deepCopy() : bb;
   }
   basicBuffer subRange(std::ptrdiff_t rangeBegin, std::ptrdiff_t rangeEnd, bool deepCopy = false){
      if((rangeBegin < 0) || rangeEnd > len_){
         return basicBuffer(nullptr, 0);
      }
      basicBuffer bb(dat_ + rangeBegin, rangeEnd - rangeBegin);
      return deepCopy ? bb.deepCopy() : bb;
   }
   T* data() { 
      return dat_; 
   } 

private:
   

   T* alloc(std::ptrdiff_t elemtCount){  
      return (elemtCount > 0 && elemtCount <= PTRDIFF_MAX) ? static_cast<T*>(allocMem(paddedLen(elemtCount * sizeof(T), cacheLineSize), cacheLineSize)) : nullptr;
   }
   //In C++17 compliant compilers, you can use std::hardware_destructive_interference_size.  
   static constexpr int cacheLineSize{64};  
   T* dat_; 
   std::ptrdiff_t len_;
   bool ownsmem_;
};

#endif

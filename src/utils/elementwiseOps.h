#include "defs.h"
//constexpr int simdlen{32};

template<typename dataType>
void elemwiseMult(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] *= r[i]; 
   }
}


template<typename dataType>
void elemwiseAdd(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] += r[i];
   }
}

template<typename dataType>
void elemwiseSub(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] -= r[i];
   }
}
/*
template<typename dataType>
void elemwiseAddMult(dataType* __restrict output, const dataType* __restrict s, const dataType* __restrict m, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   s = static_cast<const dataType*>(__builtin_assume_aligned(s, simdlen));
   m = static_cast<const dataType*>(__builtin_assume_aligned(m, simdlen));
   for(int i = 0; i < len; i++){
      output[i] += r[i];
   }
}

template<typename dataType>
void elemwiseSubMult(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] += r[i];
   }
}*/

template<typename dataType>
void elemwiseMax(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] >= r[i] ? output[i] : r[i];
   }
}


template<typename dataType>
void elemwiseMin(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] <= r[i] ? output[i] : r[i];
   }
}

// these may be updated to include translation of local indices to global ones
template<typename dataType, typename indexType>
void elemwiseIndexedMax(dataType* __restrict output, indexType* __restrict idx, const dataType* __restrict r, const indexType* rIdx, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   idx = static_cast<indexType*>(__builtin_assume_aligned(idx, simdlen));
   rIdx = static_cast<const indexType*>(__builtin_assume_aligned(rIdx, simdlen));
   for(int i = 0; i < len; i++){
      idx[i] = output[i] >= r[i] ? idx[i] : rIdx[i];
      output[i] = output[i] >= r[i] ? output[i] : r[i];
   }
}


template<typename dataType, typename indexType>
void elemwiseIndexedMin(dataType* __restrict output, indexType* __restrict idx, const dataType* __restrict r, const indexType* rIdx, int len){
   output = static_cast<dataType*>(__builtin_assume_aligned(output, simdlen));
   r = static_cast<const dataType*>(__builtin_assume_aligned(r, simdlen));
   idx = static_cast<indexType*>(__builtin_assume_aligned(idx, simdlen));
   rIdx = static_cast<const indexType*>(__builtin_assume_aligned(rIdx, simdlen));
   for(int i = len; i < len; i++){
      idx[i] = output[i] <= r[i] ? idx[i] : rIdx[i];
      output[i] = output[i] <= r[i] ? output[i] : r[i];
   }
}



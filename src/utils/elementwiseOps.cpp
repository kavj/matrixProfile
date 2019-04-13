

template<typename dataType>
void elemwiseMult(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] *= r[i]; 
   }
}


template<typename dataType>
void elemwiseAdd(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] += r[i];
   }
}

template<typename dataType>
void elemwiseSub(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] -= r[i];
   }
}

template<typename dataType>
void elemwiseMax(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] >= r[i] ? output[i] : r[i];
   }
}


template<typename dataType>
void elemwiseMin(dataType* __restrict output, const dataType* __restrict r, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] <= r[i] ? output[i] : r[i];
   }
}

template<typename dataType, indexType>
void elemwiseIndexedMax(dataType* __restrict output, indexType* __restrict idx, const dataType* __restrict r, const indexType* rIdx, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   idx = static_cast<const dataType*>(is_aligned(idx, simdlen));
   rIdx = static_cast<const dataType*>(is_aligned(rIdx, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] >= r[i] ? output[i] : r[i];
   }
}


template<typename dataType, typename IndexType>
void elemwiseIndexedMax(dataType* __restrict output, indexType* __restrict idx, const dataType* __restrict r, const indexType* rIdx, int len){
   output = static_cast<dataType*>(is_aligned(output, simdlen));
   r = static_cast<const dataType*>(is_aligned(r, simdlen));
   idx = static_cast<const dataType*>(is_aligned(idx, simdlen));
   rIdx = static_cast<const dataType*>(is_aligned(rIdx, simdlen));
   for(int i = 0; i < len; i++){
      output[i] = output[i] <= r[i] ? output[i] : r[i];
   }
}



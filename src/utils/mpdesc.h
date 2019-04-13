#include<cstddef>
#include "buffer_Iterator7.h"

template<typename dataType>
class seqArray{
   using arrayRef = basicBuffer<dataType>;
public:
   seqArray(std::ptrdiff_t seqCount, std::ptrdiff_t seqLen) : totalLen_(paddedLen(simdLen, seqLen) * seqCount), subseqs_(arrayRef(totalLen_), seqLen_(seqLen){}
   
   // simple strided iterator.
   using iterator = iterStrided<dataType>;

   iterator begin(){
      return iterator(seqs_.data(), stride_);
   }

   iterator end(){
      return iterator(seqs_.data() + totalLen_, stride_);
   }

   init(const arrayRef& ts, std::std::ptrdiff_t stride){
      // copy in as is. This will take seqCount (as defined in the constructor) sequences of seqLen_
      // Each sequence is stride_ members apart from the last as indexed with respect to ts, and these may overlap. 
      // Each sequence is copied to an aligned boundary if an aligned allocator is available. 
      // If you want to copy to a portion of this for some reason, make a subview and copy to that.

     // for(
   }

   initMeanCenter(const arrayRef& inputSeq, const arrayRef& mu, std::ptrdiff_t interSeqStep){
      // same as previous version, but also mean center each sequence, as seq - mean(seq)
   }

   initmcNormScale(const arrayRef& inputSeq, const arrayRef& mu, const arrayRef& norm, std::ptrdiff_t interSeqStep, bool isRecriprocalNorm = false){ 
      // same as above but also normalize to unit length with respect to the mean centered sequence. This is equivalent to z-normalized / sqrt(seqLen_)

   }

private:
   const std::ptrdiff_t totalLen_;
   const std::ptrdiff_t seqLen_;
   const std::ptrdiff_t stride_;
   arrayRef seqs_;
};



template<dataType>
class seqdesc{
   // maybe work in initializer List instead?
   using arrayRef = basicBuffer<dataType>;

   tsdesc(arrayRef& ts, std::ptrdiff_t windowLen = 1) : iterableLen_(ts.len() - windowLen + 1),
                                                        windowLen_(windowLen)
                                                        timeSeries_(ts), 
                                                        mu_(iterableLen_), 
                                                        inverseNorm_(iterableLen_),
                                                        df_(iterableLen_), 
                                                        dg_(iterableLen_), 
   //tsdesc(& ts, arrayRef ) : timeSeries_(ts), mu_(nullptr), inverseNorm_(nullptr), df_(nullptr), dg_(nullptr) {}
   
  
   // we need some data structure to store the buffered query entries as well? 
   // should just be able to use a basic tile with appropriate pad/striding
   std::ptrdiff_t iterLen() const{
      return iterableLen_; //timeSeries_.len() - windowLen_ + 1;
   }

   bool valid() const{
      return timeSeries_.valid() && mu_.valid() && inverseNorm.valid() && df_.valid() && dg_.valid();
   }
  
   bool init(bool needDiff_eqs = true){ 
      if(!mu_.valid() || !inverseNorm_.valid()){
         return false;
      }
      if(needDiff_eqs){
         diff_f_ = arrayRef(iterableLen_);
         diff_g_ = arrayRef(iterableLen_);  
         if(!diff_f_.valid() || !diff_g_.valid()){
            return false;
         }
      }
      init_muinvn(timeSeries_.data(), mu_.data(), inverseNorm_.data(), timeSeries_.len(), windowLen_);
      if(needDiff_eqs){
         init_dfdg(timeSeries_.data(), mu_.data(), diff_f_.data(), diff_g_.data());      
      }
   }

   seqdesc subViewLower(std::ptrdiff_t iterViewLen){
      // the arrayRefs are basicBuffers, where only the first one created owns its memory. 
      // This means that when we create a subview here, its data members (intentionally) do not own their underlying buffer space
      // This allows us to have multiple readable views without extensive copying. It might be a good idea to make a const iterator for this.
      if(iterViewLen <= iterableLen_){
         return seqdesc(*this, 0, iterViewLen_);
      }
      return seqdesc(0, 0); // null 
   }
   
   seqdesc subViewUpper(std::ptrdiff_t iterableLen){
      if(iterableLen <= iterableLen_ && valid()){
         return seqdesc(*this, iterableLen_ - iterableLen, iterableLen_); 
      }
      return seqdesc(0, 0);
   }

   seqArray meanCenteredSet(){
      
   }

   seqArray meanCenteredNormSet(){
      
   }

private:
   std::ptrdiff_t iterableLen_; // by default, we attenuate edges and leaves this as len(timeSeries_) - windowLen_ + 1
   std::ptrdiff_t windowLen_;   // might add non-attenuated options later
   arrayRef timeSeries_;
   arrayRef mu_;
   arrayRef inverseNorm_;
   arrayRef diff_f_;  // difference equations, used to update covariance
   arrayRef diff_g_;
};


template<typename dataType, typename indexType>
class indexedProfile{
   using arrayRef = basicBuffer<dataType>;
   using indexRef = basicBuffer<indexType>;
public:
   // convenience functions to merge 2 profiles. offsetIdxBy means that we need to apply
   // an offset to the incoming index members prior to merging them. This should call a simd friendly function. 
   void elemtwiseMax(const indexedProfile& idxProf, std::ptrdiff_t offsetIdxBy = 0){
      
   }

   void elemtwiseMin(const indexedProfile& idxProf, std::ptrdiff_t offsetIdxBy = 0){

   }
   
private:
   arrayRef p_;
   indexRef pIndex_;
};


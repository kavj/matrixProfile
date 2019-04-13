#ifndef ITERS
#define ITERS

//#include<cstddef>

// Since padding introduces its own problems, strided iterator store explicit iterable length to avoid undefined behavior
// If you need an out of place copy of the underlying buffer, use a dereference followed by a deep copy.

template<typename T>
class iterStrided{  
public:
   iterStrided(T* dat, std::ptrdiff_t iterableLen, std::ptrdiff_t stride) : dat_(dat), iterableLen_(iterableLen), stride_(stride){}
   std::ptrdiff_t operator - (const iterStrided<T>& t) const{ 
      return dat_ - t.dat_;
   }
   bool operator !=(const iterStrided<T>& t) const{
      return dat_ != t.dat_; 
   }
   bool operator ==(const iterStrided<T>& t) const{
      return dat_ == t.dat_;
   } 
   std::ptrdiff_t iterableLen() const{
      return iterableLen_;
   }
   iterStrided<T> operator++(int) = delete;
   iterStrided<T> operator++(){
      if(iterableLen_ == 0){ // If we are already 1 past the last element, meaning "end", return nullptr. 
         return iterStrided<T>(nullptr, stride_);
      } 
      std::ptrdiff_t step = stride_ <= iterableLen_  ? stride_ : iterableLen_;
      dat_ += step;
      iterableLen_ -= step; 
      return *this;
   }
   T* operator*(){
      return iterableLen_ > 0 ? dat_ : nullptr; 
   }

private:
   T* dat_;
   std::ptrdiff_t iterableLen_;
   const std::ptrdiff_t stride_;
};


#endif

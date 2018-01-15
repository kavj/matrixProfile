#include "../arch/avx2arith_2.hpp"

// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al

// The first argument is always updated with the resulting sum or product and the second is updated to reflect the error term.
using namespace vmth;
#define fact 2 << 27 


// not sure if it's legal to use an enumerated type in second arg
template<typename dtype>
static inline double factor(dtype a){
   // this needs to be made compliant with vector types, perhaps by overloading set(scalar) for scalar operands
   //dtype b = static_cast<dtype>(a);
   static const dtype c = bcast(static_cast<dtype>(fact));
   return a*c;
//   return _mm256_castsi256_pd(bcast(fact))*a;
}

template<typename dtype>
static inline void xadd(dtype &a, dtype &b){
   dtype c = a + b;
   dtype d = c - a;
   b = ((a - (c - d)) + (b - d));
   a = c;
}

template<typename dtype>
static inline void xsplit(dtype &a,dtype &b){
   dtype f = setzero(); 
   
//   __m256d f = factor<__m256d>(a);
   b = (f - (f - a));
   a = a - f;
}

// this should probably move to a separate header or the top or something
template<typename dtype>
static inline void xmul(dtype &a,dtype &b){
#ifdef __FMA__
   dtype c = a*b;
   b = fms(a,b,c);
   a = c;
#else
   dtype c = a;
   a *= b;
   dtype f = b;
   dtype d;
   dtype g;
   xsplit(c,d);
   xsplit(f,g);
   b = d*g - (((a - c*f) - c*f) - c*g);
#endif
}


template<typename dtype>
void xsum(dtype *a, dtype &s, dtype &e, int len){
   for (int i = 0; i < len; i++){
      dtype b = a[i];
      xadd(s,b);
      e += b;
   }
}

template<typename dtype, typename vtype>
static inline vtype xsInv_kern(dtype *a, dtype *mu, dtype *sI, int offset, long winlen){
   vtype z = uload(mu,offset);
   vtype p = uload(a,offset) - z;
   vtype s = p;
   xmul(p,s);
   dtype u = static_cast<dtype>(winlen);
   vtype t = bcast(u);
 //  ustore(t,sI,0);
   for(int j = offset+1; j < offset+winlen; j++){
       vtype h = uload(a,j) - z;
       vtype r = h;
       xmul(h,r);
       xadd(p,h);
       s += (h+r);
   }
   p += s;
  
   return p;
}


template<typename dtype>
void xsInv(dtype *a, dtype *mu, dtype *sI, int len, int winlen){
   
   int k = len - winlen + 1;
   long z = static_cast<long>(winlen);
   k -= k%8;
   #pragma omp parallel for
   for(int i = 0; i < k; i+=8){
      ustore(xsInv_kern<dtype,__m256d>(a,mu,sI,i,z),sI,i);
      ustore(xsInv_kern<dtype,__m256d>(a,mu,sI,i+4,z),sI,i+4);
   }
//   if(k < (len-winlen+1){

//   }
}


template<typename dtype>
void xmean_windowed(dtype *a, dtype *output, int len, int winlen){
   int mlen = len - winlen + 1;
   dtype s = 0;
   dtype e = 0;
   xsum(a,s,e,winlen);
   output[0] = (s + e)/winlen;
   for (int i = winlen; i < len; i++){
      dtype b = -1*a[i-winlen];
      xadd(s,b);
     // xsub(s,b);
      e += b; 
      b = a[i];
      xadd(s,b);
      e += b;
      output[i-winlen+1] = (s + e)/winlen;
   }
}



#include "../mp/avx2arith_2.hpp"

// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al

// The first argument is always updated with the resulting sum or product and the second is updated to reflect the error term.

#define fact 134217728 


// not sure if it's legal to use an enumerated type in second arg
template<typename dtype>
static inline double factor(dtype a){
   // this needs to be made compliant with vector types, perhaps by overloading set(scalar) for scalar operands
   return fact*a;
}

template<typename dtype>
static inline void xadd(dtype &a, dtype &b){
   dtype c = a + b;
   dtype d = c - a;
   b = (a - (c - d)) + (b - d);
   a = c;
}

template<typename dtype>
static inline void xsplit(dtype &a,dtype &b){
   dtype f = factor(a);
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
static inline dtype xadd(dtype &a, dtype &b){
   dtype c = a + b;
   dtype d = c - a;
   b = (a - (c - d)) + (b - d);
   a = c;
}

// double check this one
template<typename dtype>
static inline dtype xsub(dtype &a, dtype &b){
   dtype c = a - b;
   dtype d = c - a;
   b = (a - (c - d)) - (b + d);
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
void xssq(dtype *a, dtype *b, dtype *output, int winlen){
   int mlen = n - winlen + 1;
   vtype z = bcast(b,0);
   vtype p = uload(a,0) - z;
   vtype s = x;
   xmul(x,y);
   for(int i = 1; i < winlen; i++){
      vtype h = uload(a,i) - z;
      vtype r = h;
      xadd(x,
   }
}


// simd isn't worthwhile with a simple windowed sum method, particularly because it makes error analysis harder
template<typename dtype>
void xmean_windowed(dtype *a, dtype *output, int len, int winlen){
   dtype s = 0;
   dtype e = 0;
   xsum(a,s,e);
   vtype = bcast(winlen);
   output[0] = (s + e)/winlen;
   for (int i = 1; i < len-winlen+1; i++){
      dtype b = a[i];
      xsub(s,b);
      e += b; 
      b = a[i+winlen];
      xadd(s,b);
      e += b;
      output[i] = (s + e)/winlen;
   }
}


// this should mainly  be vectorized if we're allowed to use fused multiply add
template<typename dtype, typename vtype>
void xvar_windowed(dtype *a, dtype *b, int len, int wlen){
   vtype s = setzero(); // setzero();
   vtype e = setzero();   //setzero();
   for(int i = 0; i < wlen; i++){
      
   }

}



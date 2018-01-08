

// Collection of in place 2 operand arithmetic and sliding window sum and average. These are based on 
// ACCURATE SUM AND DOT PRODUCT Ogita et al

// The first argument is always updated with the resulting sum or product and the second is updated to reflect the error term.


// not sure if it's legal to use an enumerated type in second arg
template<typename dtype>
static inline double factor(dtype a){
   // this needs to be made compliant with vector types, perhaps by overloading set(scalar) for scalar operands
   return set(factor)*a;
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
#ifdef __FMA__
template<typename dtype>
static inline void xmul(dtype &a, dtype &b){
   dtype a_1 = a;
   dtype a_2;
   a = a*b;
   dtype b_1 = b;
   dtype b_2;
   xsplit(a_1,a_2);
   xsplit(b_1,b_2);
   b = a_2*b_2 - (((a - a_1*b_1) - a_2*b_1) - a_1*b_2);
}
#else
template<typename dtype>
static inline void xmul(dtype &a,dtype &b){
   dtype c = a*b;
   b = fms(a,b,c);
   a = c;
}
#endif

template<typename dtype>
static inline dtype xadd(dtype &a, dtype &b){
   dtype c = a + b;
   dtype d = c - a;
   b = (a-(c-d))+(b-d);
   a = c;
}

// double check this one
template<typename dtype>
static inline dtype xsub(dtype &a, dtype &b){
   dtype c = a - b;
   dtype d = c - a;
   b = (a - (c - d)) - (b + d);
}

template<typename dtype,typename vtype>
dtype xsum(dtype *a, int len){
   dtype b = 
   
}

// vtype is vector type, indicating what type of vector to use
// this makes no alignment assumptions on the input data, only output
template<typename dtype, typename vtype, int stride>
void xsum_windowed(dtype *a, dtype *output, int len, int winlen){
   vtype s = uload(a,0);
   vtype e = set(0); 
   for (int i = stride; i < winlen; i+= stride){
      vtype b = uload(a,i);
      xadd(s,b);
      e = e + b;
   }
   astore(s + e,output,0);
   for (int i = winlen; i < len; i += stride){
      vtype b = uload(a,i-winlen);
      
   }
}




#ifndef MAX_REDUCE
#define MAX_REDUCE
#include "../arch/avx256.h"
using namespace avx256_t;


struct rpair { __m256d val; __m256i index;};

 
static inline struct rpair __attribute__((always_inline)) max_reduce_8x1(const __m256d& r0, const __m256d& r1, const __m256d& r2, const __m256d& r3,  const __m256d& r4, const __m256d& r5, const __m256d& r6, const __m256d& r7){

   __m256i mask1, mask2, mask3, mask4, mask5, mask6;
   __m256d  op1, op2, op3, op4;

   mask1 = r0 > r1;
   mask2 = r2 > r3;
   mask3 = r4 > r5;
   mask4 = r6 > r7;

   const __m256i four = brdcst(4);
   const __m256i eight = brdcst(8);
  
   op1 = vmax(r0,r1);
   op2 = vmax(r2,r3);
   op3 = vmax(r4,r5);
   op4 = vmax(r6,r7);
  
   __m256i g = set(0,1,2,3); 
   __m256i h = g + four;
   mask1 = blend(g, h, mask1);
   g += eight;
   h += eight;
   mask2 = blend(g, h, mask2);
   g += eight;
   h += eight;
   mask3 = blend(g, h, mask3);
   g += eight;
   h += eight;
   mask4 = blend(g, h, mask4);
   
   mask5 = op1 > op2;
   mask6 = op3 > op4;
   
   op1 = vmax(op1, op2);
   op3 = vmax(op3, op4);

   mask1 = blend(mask1, mask2, mask5);
   mask3 = blend(mask3, mask4, mask6);
   
   struct rpair x;
   x.val = vmax(op1, op3);
   x.index = blend(mask1, mask3, op1 > op3);
   return x;
}

 
static inline struct rpair __attribute__((always_inline)) max_pair(const struct rpair& r0, const struct rpair& r1){
   struct rpair x;
   x.index = blend(r0.index, r1.index, (__m256i)(r0.index < r1.index));
   x.val = blend(r0.val, r1.val, r0.val < r1.val);
   return x;
}



#endif

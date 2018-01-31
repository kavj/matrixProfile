





void mpkern(double* __restrict__ cx, const double*  __restrict__ dx, const double* __restrict__ xf, const double* __restrict__ s, double* __restrict__ buf, int diag, int offset); //__attribute__ ((always_inline));



void mpkern2(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset);

static inline void reduce(double* __restrict__ buf, double*  __restrict__ output, double* __restrict__ outputi, int diag, int k) __attribute__ ((always_inline));



/*
void mpkern(double* __restrict__ cx, const double* __restrict__ dx, const double* __restrict__ xf, double* __restrict__ s, double* __restrict__ buf, int diag, int offset){
   cx = (double*)__builtin_assume_aligned(cx,64);
   dx = (double*)__builtin_assume_aligned(dx,64);
   xf = (double*)__builtin_assume_aligned(xf,64);
   s  = (double*)__builtin_assume_aligned(s,64);
   buf =(double*)__builtin_assume_aligned(buf,64);

   int m = diag+offset;
   for(int i = 0; i < 1024; i++){
      double c = cx[i];
      double s1 = s[i];
      c = c + dx[i]*xf[i+diag];
      c = c + dx[i+diag]*xf[i];
      double s2 = s[i];
      double d = c*s1*s2;
      buf[i+64] = d;
   }
}
*/

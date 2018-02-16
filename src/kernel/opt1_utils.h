

template<typename dt>
static inline void fma_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt &ra, dt *a, int index) __attribute__ ((always_inline));

template<typename dt>
static inline void prod_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt &ra, dt *a, int index) __attribute__ ((always_inline));

template<typename dt>
static inline void load_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *a, int index) __attribute__ ((always_inline));

template<typename dt>
static inline void load_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *a, int index) __attribute__ ((always_inline));

//static inline void fma_x8<dt>(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt &ra, dt *a, int index) __attribute__ ((always_inline));


template<typename dt>
static inline void load_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *a, int index){
   ca = a[index];
   cb = a[index+1];
   cc = a[index+2];
   cd = a[index+3];
   ce = a[index+4];
   cf = a[index+5];
   cg = a[index+6];
   ch = a[index+7];
}

template<typename dt>
static inline void store_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *a, int index){
   a[index] = ca;
   a[index+1] = cb;
   a[index+2] = cc;
   a[index+3] = cd;
   a[index+4] = ce;
   a[index+5] = cf;
   a[index+6] = cg;
   a[index+7] = ch;
}

template<typename dt>
static inline void fma_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt &ra, dt *a, int index){
   ca = fma(ca,ra,a[index]);
   cb = fma(cb,ra,a[index+1]);
   cc = fma(cc,ra,a[index+2]);
   cd = fma(cd,ra,a[index+3]);
   ce = fma(ce,ra,a[index+4]);
   cf = fma(cf,ra,a[index+5]);
   cg = fma(cg,ra,a[index+6]);
   ch = fma(ch,ra,a[index+7]);
}


template<typename dt>
static inline void prod_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *a, int index){
   ca *= a[index];
   cb *= a[index+1];
   cc *= a[index+2];
   cd *= a[index+3];
   ce *= a[index+4];
   cf *= a[index+5];
   cg *= a[index+6];
   ch *= a[index+7];
}


template<typename dt>
static inline void prod_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt &ra){
   ca *= ra;
   cb *= ra;
   cc *= ra;
   cd *= ra;
   ce *= ra;
   cf *= ra;
   cg *= ra;
   ch *= ra;
}

template<typename dt>
static inline void prod_oop_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch,
                                   dt &da, dt &db, dt &dc, dt &dd, dt &de, dt &df, dt &dg, dt &dh, dt* q, int index){
   da = ca * q[index];
   db = cb * q[index+1];
   dc = cc * q[index+2];
   dd = cd * q[index+3];
   de = ce * q[index+4];
   df = cf * q[index+5];
   dg = cg * q[index+6];
   dh = ch * q[index+7];
}

template<typename dt, typename dti>
static inline void min_x8(dt &ca, dt &cb, dt &cc, dt &cd, dt &ce, dt &cf, dt &cg, dt &ch, dt *p, dti *pi, int index, int offset){
   if(ca < p[index]){
      p[index] = ca;
      pi[index] = offset;
   }

}

template<typename dt>
static inline void collapse_x8(){

}




 

namespace scalar{
#ifdef __GNUC__
#define wrapper __attribute__((always_inline,artificial))
#else
#define wrapper
#endif



#define stride 1

static inline double wrapper aload(const double *a, int offset){
   return a[offset];
}

static inline __m256d wrapper uload(const double *a, int offset){
   return a[offset];
}

static inline void wrapper astore(const double& b,double *a, int offset){
   a[offset] = b;
}

static inline void wrapper ustore(const double& b, double* a,int offset){
   a[offset] = b;
}

static inline double wrapper mul_add(const double &a,const double &b,const double &c){
   return  a*b + c;
}

static inline double wrapper mul_sub(const double &a, const double &b,const double &c){
   return a*b - c;
}

static inline double wrapper mul_nadd(const double &a, const double &b, const double &c){
   return c - a*b;
}

static inline double wrapper vmax(const double &a, const double &b){
   return max(a,b);
}

static inline double wrapper vmin(const double &a, const double &b){
   return min(a,b);
}

static inline long wrapper aload(const long* a, int offset){
   return  a[offset];
}

static inline int wrapper aload(const int* a, int offset){
   return a[offset];
}

static inline long wrapper uload(const long* a, int offset){
   return a[offset];
}

static inline int wrapper uload(const int* a, int offset){
   return a[offset];
}

static inline void wrapper astore(const long &oper, long* a, int offset){
   a[offset] = oper;
}

static inline void wrapper astore(const int &oper, int* a, int offset){
   a[offset] = oper;
}

static inline void wrapper ustore(const long &oper, long* a,int offset){
   a[offset] = oper;
}

static inline void wrapper ustore(const int &oper, int* a,int offset){
   a[offset] = oper;
}

static inline double wrapper brdcst(const double* a,int offset){
   return a[offset];
}

static inline long wrapper brdcst(const long a){
   return a;
}

static inline int wrapper brdcst(const int a){
   return a;
}

static inline __m256d wrapper brdcst(const double a){
   return a;
}
}


#endif
#endif

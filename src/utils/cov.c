#define prefalign 64 
// Todo: Test whether simplified indexing would break compiler optimizations
#define wid 64 

void center_query(const double* restrict ts, const double* restrict mu, double* restrict q, int sublen){
   q = (double*)__builtin_assume_aligned(q,prefalign);
   const int aligned = sublen <= wid ? sublen : sublen - sublen % wid;
   for(int i = 0; i < aligned; i+= wid){
      for(int j = i; j < i + wid; j++){
         q[j] = ts[j] - mu[0];
      }
   }
   for(int i = aligned; i < sublen; i++){
      q[i] = ts[i] - mu[0];
   }
}

void batchcov(const double* restrict ts, const double* restrict mu, const double* restrict query, double* restrict cov, int count, int sublen){
   query = (double*)__builtin_assume_aligned(query,prefalign);
   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   const int alcount = count <= wid ? count : count - count % wid;
   for(int i = 0; i < alcount; i+= wid){
      for(int j = 0; j < wid; j++){
         cov[i + j] = 0;
      }
      for(int k = 0; k < sublen; k++){
         for(int j = 0; j < wid; j++){
            cov[i + j] += (ts[i + j + k] - mu[i + j]) * query[k];
         }
      }
   }
   for(int i = alcount; i < count; i++){
      cov[i] = 0;
   }
   for(int i = 0; i < sublen; i++){
      for(int j = alcount; j < count; j++){
         cov[j] += (ts[i + j] - mu[j]) * query[i]; 
      }
   }
}

#define prefalign 64 
// Todo: Test whether simplified indexing would break compiler optimizations
#define wid 64 

void center_query(const double* __restrict__ ts, const double* __restrict__ mu, double* __restrict__ q, int sublen){
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

void batchcov(const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ query, double* __restrict__ cov, int count, int sublen){
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

// Will update these later, They will eventually displace xprec, which is still a bit sloppy 
// Todo: verify that zero length windows are accounted for prior to reaching this point
#include<array>

void trailing_mean(double* __restrict a, double* __restrict mu, long long len, long long winlen){
   double accum = a[0];
   double resid = 0;
   for(long long i = 1; i < winlen; i++){
      double m = a[i];
      double p = accum;
      accum += m;
      double q = accum - p;
      resid += ((p - (accum - q)) + (m - q));
   }
   mu[0] = accum + resid;
   for(long long i = winlen; i < len; i++){
      double m = a[i - winlen];
      double n = a[i];
      double p = accum - m;
      double q = p - accum;
      double r = resid + ((accum - (p - q)) - (m + q));
      accum = p + n;
      double s = accum - p;
      resid = r + ((p - (accum - s)) + (n - s));
      mu[i - winlen + 1] = accum + resid;
   }
   for(long long i = 0; i < len - winlen; i++){
      mu[i] /= winlen;
   }
}

void trailing_inverse_centered_norm(double* __restrict a, double* __restrict mu, double* __restrict invn, long long len, long long winlen){
   const long long alignedwinlen = winlen - winlen % 128;
   #pragma omp parallel for
   for(long long i = 0; i < len - winlen; i++){
      std::array<double, 128> resid;
      std::array<double, 64> aux;
      for(long long j = i; j < i + alignedwinlen; j += 64){
         for(long long k = j; k < j + 64; k++){
            aux[k - j] = a[j] - mu[i];
     	    //aux[k - j] *= aux[k - j];
	 }
	 for(long long k = j; k < j + 64; k++){
           invn[k] = aux[k - j] * aux[k - j]; 
	 }
	 for(long long k = j; k < j + 128; k++){

	 }
      }
      for(long long j = i + alignedwinlen; j < winlen; j++){
         for(long long k = j; k < j + 64; k++){
            aux[k - j] = a[j] - mu[i];
     	    //aux[k - j] *= aux[k - j];
	 }
	 for(long long k = j; k < j + 64; k++){
           invn[k] = aux[k - j] * aux[k - j]; 
	 }
	 for(long long k = j; k < j + 128; k++){

	 }
      }
      for(long long j = i + alignedwinlen; j < winlen; j++){
        // nosimd 
      }
   }
}


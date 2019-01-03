#include<cmath>
#include<array>
constexpr int prefalign  = 64;
constexpr int blocklen = 32;


void center_query(const double* __restrict ts, const double* __restrict mu, double* __restrict query, int winlen){
   query = (double*)__builtin_assume_aligned(query, prefalign);
   mu = (double*) __builtin_assume_aligned(mu, prefalign);
   #pragma omp simd aligned(query, mu : prefalign) safelen(128)  
   for(int i = 0; i < winlen; i++){
      query[i] = ts[i] - mu[0];
   }
}


#pragma omp declare simd aligned(ts, mu, query, cov : prefalign) 
void batchcov(const double* __restrict ts, const double* __restrict mu, const double* __restrict query, double* __restrict cov, int count, int winlen){
   ts = (double*) __builtin_assume_aligned(ts,prefalign);
   mu = (double*) __builtin_assume_aligned(mu,prefalign);
   query = (double*)__builtin_assume_aligned(query,prefalign);
   cov = (double*)__builtin_assume_aligned(cov,prefalign);
   for(int i = 0; i < count; i++){
      cov[i] = 0;
      #pragma omp simd aligned(ts, mu, query, cov : prefalign) safelen(128)
      for(int j = 0; j < winlen; j++){
         cov[i + j] += (ts[i + j] - mu[i]) * query[j];
      }
   }
}




// The following 2 functions are based on the work in Ogita et al, Accurate Sum and Dot Product

void sw_mean(double* __restrict a, double* __restrict mu, int len, int winlen){
   double accum = a[0];
   double resid = 0;
   for(int i = 1; i < winlen; i++){
      double m = a[i];
      double p = accum;
      accum += m;
      double q = accum - p;
      resid += ((p - (accum - q)) + (m - q));
   }
   mu[0] = accum + resid;
   for(int i = winlen; i < len; i++){
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
   for(int i = 0; i < len - winlen + 1; i++){
      mu[i] /= winlen;
   }
}


void sw_inv_meancentered_norm(double* __restrict a, double* __restrict mu, double* __restrict invn, int len, int winlen){
   mu = static_cast<double*>(__builtin_assume_aligned(mu, prefalign));
   invn = static_cast<double*>(__builtin_assume_aligned(invn, prefalign));
   const int alignedwinlen = (winlen >= blocklen) ? (winlen - winlen % blocklen) : 0;
   #pragma omp parallel for
   for(int i = 0; i < len - winlen + 1; i++){
      double accsum = 0;
      double accresid = 0;
      for(int j = i; j < i + alignedwinlen; j += blocklen){
         std::array<double, blocklen> sq;
	 std::array<double, blocklen> aux;
         for(int k = 0; k < blocklen; k++){
            aux[k] = a[j + k] - mu[i];
	 }
	 for(int k = 0; k < blocklen; k++){
            sq[k] = aux[k] * aux[k];
	 }
	 for(int k = 0; k < blocklen; k++){
            aux[k] = fma(aux[k], aux[k], -sq[k]);
	 }
	 for(int k = 0; k <  blocklen; k++){
     	    double term  = sq[k]; 
	    double priorsum = accsum;
            accsum += term;
            double r = accsum - priorsum;
	    double resid = (priorsum - (accsum - r) + (term - r));
	    accresid += (aux[k] + resid);
	 }
      }
      for(int j = i + alignedwinlen; j < i + winlen; j++){
         double mc = a[j] - mu[i];
	 double term = mc * mc;
	 double residsq = fma(mc, mc, -term);
         double priorsum = accsum;  
         accsum += term;
	 double r = accsum - priorsum;
	 double residadd = (priorsum - (accsum - r) + (term - r));
	 accresid += (residsq + residadd);
      }
      invn[i] = accsum + accresid;
   }
   // debating whether to split this section
   #pragma omp parallel for
   for(int i = 0; i < len - winlen + 1; i++){
      invn[i] = 1/sqrt(invn[i]);
   }
}


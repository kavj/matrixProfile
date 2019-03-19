#include<cmath>
#include<array>
#include "cov.h"
constexpr int simdAlignment  = 64;
constexpr int blocklen = 32;
constexpr int wid = 128; // optimize out later
constexpr int prefalign = 64; // redundant, get rid of later

void centerQuery(const double* __restrict ts, double* __restrict qBuffer, const double mu, int windowLen){
   qBuffer = (double*)__builtin_assume_aligned(qBuffer, simdAlignment);
   #pragma omp simd aligned(qBuffer : simdAlignment)
   for(int i = 0; i < windowLen; i++){
      qBuffer[i] = ts[i] - mu;
   }
}


void crossCov(const double* __restrict__ ts, const double* __restrict__ mu, const double* __restrict__ query, double* __restrict__ cov, int count, int windowLen){
   query = (double*)__builtin_assume_aligned(query,simdAlignment);
   cov = (double*)__builtin_assume_aligned(cov,simdAlignment);
   const int alcount = count <= wid ? count : count - count % wid;
   for(int i = 0; i < alcount; i+= wid){
      for(int j = 0; j < wid; j++){
         cov[i + j] = (ts[i + j] - mu[i + j]) * query[0];
      }
      for(int k = 1; k < windowLen; k++){
         for(int j = 0; j < wid; j++){
            cov[i + j] += (ts[i + j + k] - mu[i + j]) * query[k];
         }
      }
   }
   for(int i = alcount; i < count; i++){
      cov[i] = (ts[i] - mu[i]) * query[0];
   }
   for(int i = 1; i < windowLen; i++){
      for(int j = alcount; j < count; j++){
         cov[j] += (ts[i + j] - mu[j]) * query[i]; 
      }
   }
}


void windowedMean(double* __restrict a, double* __restrict mu, int len, int winlen){
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


void windowedInverseCenteredNorm(double* __restrict a, double* __restrict mu, double* __restrict invn, int len, int winlen){
   mu = static_cast<double*>(__builtin_assume_aligned(mu, prefalign));
   invn = static_cast<double*>(__builtin_assume_aligned(invn, prefalign));
   const int alignedwinlen = (winlen >= blocklen) ? (winlen - winlen % blocklen) : 0;
  // #pragma omp parallel for
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
   //#pragma omp parallel for
   for(int i = 0; i < len - winlen + 1; i++){
      invn[i] = 1/sqrt(invn[i]);
   }
}


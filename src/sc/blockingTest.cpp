#define _POSIX_C_SOURCE 200112L 
//#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <string>
//#include <array>
#include <immintrin.h>
#include <assert.h>



/*void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output,int n,int m){
    double a[128];
    unsigned long b = 0;
    //printf("%d  %d \n",n,m);
    for(int i = 0; i < (n-m); i +=64){
        for(int j = i; j < n-m;j+=4){
            __m256d c = _mm256_load_pd(Cxy+j);
            for(int k = 0; k < 64; k+=4){
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_loadu_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_loadu_pd(dX+j+k+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+k,c);
                _mm256_store_pd(a+k,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k))); 
                _mm256_store_pd(a+k,_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k+m)));
                __m256d f = _mm256_load_pd(output+k);
                __m256d g = _mm256_load_pd(a+k);
                __m256d h = _mm256_loadu_pd(a+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,13);
                __m256d d2 = _mm256_cmp_pd(g,h,13);         
                g = _mm256_blendv_pd(g,f,d1);
                h = _mm256_blendv_pd(h,f,d2);
                _mm256_store_pd(output+k,g);
                _mm256_store_pd(output+k+m,h);
                b ++;
            }
        }
    }
    printf("%lu %d %d \n",b,n,n-m+1);
}*/


void  accumTest11(double* Cxy, double* dX, double* dF, double*s, double* output,int* outputI,int n,int m,int iters){
    double fwbuf[256];
    double bwbuf[256];
    unsigned long blah = 0;
    m /= 4;
    m *= 4;
    for(int iters = 0; iters < n-m; iters++){
    for(int i = 2*m; i < n-4*m; i+= m){
        __m256d fw = _mm256_load_pd(Cxy+i);
        __m256d bw = _mm256_load_pd(Cxy+i);
        for(int j = 0; j < m; j+=4){
            fw = _mm256_fmadd_pd(_mm256_load_pd(dF+i+j),_mm256_load_pd(dX+i+j),fw);
            _mm256_storeu_pd(fwbuf+4*j,fw);
        }
        for(int j = m; j > 0; j-= 4){
            bw = _mm256_fmadd_pd(_mm256_load_pd(dF+i+j),_mm256_load_pd(dX+i+j),bw);
            _mm256_store_pd(bwbuf+4*j,bw);
        }
        for(int j = 0; j < m; j+= 4){
            _mm256_store_pd(output+i+j,_mm256_add_pd(_mm256_load_pd(fwbuf+j),_mm256_load_pd(bwbuf+j)));
        }
        blah += 2*m;
    }
    }
    printf("%lu\n",blah);
}

void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output,int* outputI,int n,int m,int iters){
    //std::array<__m256d,64> a;
    double a[64];
    unsigned long t3 = 0;
    for(int p = 0; p < iters; p++){
    for(int i = 0; i < (n-m); i +=32){
        __m256d c = _mm256_load_pd(Cxy+i);
        for(int j = 0; j < i;j++){
            for(int k = 0; k < 32; k+=4){
                __m256i ind3 = _mm256_set1_epi64x((long long)i);
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_load_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_load_pd(dX+j+k+m);
                __m256i ind1 = _mm256_load_si256((__m256i*)outputI);
                __m256i ind2 = _mm256_load_si256((__m256i*)outputI+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+k,c);
                __m256d g = _mm256_mul_pd(_mm256_load_pd(s+j+k+m),_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k))); 
                __m256d f = _mm256_load_pd(output+k);
                __m256i t = _mm256_load_si256((__m256i*)outputI+j);
                __m256i t2 = _mm256_load_si256((__m256i*)outputI+j+m);
                __m256d h = _mm256_load_pd(output+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,_CMP_GT_OQ);
                __m256d d2 = _mm256_cmp_pd(g,h,_CMP_GT_OQ);         
                g = _mm256_blendv_pd(g,f,d1);
                h = _mm256_blendv_pd(g,h,d2);
                ind1 = _mm256_blendv_epi8(ind3,ind1,(__m256i)d1);
                ind2 = _mm256_blendv_epi8(ind3,ind2,(__m256i)d2);
                _mm256_store_pd(output+k,g);
                _mm256_store_pd(output+k+m,h);
                t3+=4;
            }
        }
        _mm256_store_pd(Cxy+i,c);
    }
    }
    printf("%lu \n",t3,n,n-m+1);
}


//void  accumTest4(T, dX, dF, sigmaInv, buffer, bufferI,1024*i, m,1024, m);
/*void  accumTest4(double* Cxy, double* dX, double* dF, double*s, double* output,int* outputI,int base, int offs,int n,int m){
    double a[64];
    unsigned long t3 = 0;
    for(int i = base; i < base+n-m; i += 512){
        __m256d c = _mm256_load_pd(Cxy+i);  // i indicates the lag factor in the case of Cxy (unfortunately opposite my usual conventions)
        for(int k = i; k < k+512;k+=4){     // j indicates how far forward
            for(int j = offs+k; j < offs+k+512; j++){
            for(int l = 0; l < 32; l+=4){
                __m256i ind3 = _mm256_set1_epi64x((long long)i);
                __m256d op1 = _mm256_broadcast_sd(dF+i+k);
                __m256d op2 = _mm256_load_pd(dX+i+k);
                __m256d op3 = _mm256_broadcast_sd(dF+j+k+m);
                __m256d op4 = _mm256_load_pd(dX+j+k+m);
                __m256i ind1 = _mm256_load_si256((__m256i*)outputI);
                __m256i ind2 = _mm256_load_si256((__m256i*)outputI+m);
                c = _mm256_fmadd_pd(op1,op2,c);
                c = _mm256_fmadd_pd(op3,op4,c);
                _mm256_store_pd(a+k,c);
                __m256d g = _mm256_mul_pd(_mm256_load_pd(s+j+k+m),_mm256_mul_pd(_mm256_load_pd(a+k),_mm256_loadu_pd(s+j+k))); 
                __m256d f = _mm256_load_pd(output+k);
                __m256i t = _mm256_load_si256((__m256i*)outputI+j);
                __m256i t2 = _mm256_load_si256((__m256i*)outputI+j+m);
                __m256d h = _mm256_load_pd(output+k+m);
                __m256d d1 = _mm256_cmp_pd(g,f,_CMP_GT_OQ);
                __m256d d2 = _mm256_cmp_pd(g,h,_CMP_GT_OQ);         
                g = _mm256_blendv_pd(g,f,d1);
                h = _mm256_blendv_pd(g,h,d2);
                ind1 = _mm256_blendv_epi8(ind3,ind1,(__m256i)d1);
                ind2 = _mm256_blendv_epi8(ind3,ind2,(__m256i)d2);
                _mm256_store_pd(output+k,g);
                _mm256_store_pd(output+k+m,h);
                t3+=4;
            }
        }
        _mm256_store_pd(Cxy+i,c);
      }
    printf("%lu \n",t3,n,n-m+1);
    }
} */
void scalTest4(double* Cxy, double* dX, double* dF, double*s, double* output,int n,int m){
double a;
int b = 0;
for(int i = 0; i < (n-m); i +=64){
    for(int j = i; j < n-m;j+=4){
        double c = Cxy[j];
        for(int k = 0; k < 64; k+=4){
            double op1 = dF[i+k];
            double op2 = dX[i+k];
            double op3 = dF[j+k+m];
            double op4 = dX[j+k+m];
            c += op1*op2 + op3*op4;
            a = c*s[j+k]*s[j+k+m];
            double z = output[k];
            double g = a > z ? g : z;
            output[k] = g;
            b ++;
        }
    }
}
printf("%lu %d %d \n",b,n,n-m+1);
}
   
 

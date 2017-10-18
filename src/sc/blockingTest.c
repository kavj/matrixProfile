#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <immintrin.h>




void regBlockTest(double* T, double* mu, double* output, const int minLag, int n, int m){
    double a = 0;
    double b = 0;
    double c = 0;
    double d = 0;
 /*   double e = 0;
    double f = 0;
    double g = 0;
    double h = 0;
*/ 
    //int j = 0;
    for(int j = 0; j < m; j++){
        for(int i = j; i < j+m; i++){
            
            double p0 = T[i] - mu[i];
            double p1 = T[i+1] - mu[i+1];
            double p2 = T[i+2] - mu[i+2];
            double p3 = T[i+3] - mu[i+3];
            double t0 = T[i+minLag+1]-mu[i+minLag+1];
            double t1 = T[i+minLag+2]-mu[i+minLag+2];
            double t2 = T[i+minLag+3]-mu[i+minLag+3];
            double t3 = T[i+minLag+4]-mu[i+minLag+4];
            a += p0*t0;
            b += p0*t1;
            c += p0*t2;
            d += p0*t3;
            a += p1*t1;
            b += p1*t2;
            c += p1*t3;
            //d += p1*t4;
        }
    
    //for(int j = 0; j < 2; j++){
/*    for(int i = j; i < j+m; i++){
        a += (T[i]-mu[i])*(T[i+minLag]-mu[i+minLag]);
        b += (T[i]-mu[i])*(T[i+minLag+1]-mu[i+minLag+1]);
        c += (T[i]-mu[i])*(T[i+minLag+2]-mu[i+minLag+2]);
        d += (T[i]-mu[i])*(T[i+minLag+3]-mu[i+minLag+3]);
        e += (T[i]-mu[i])*(T[i+minLag+4]-mu[i+minLag+4]);
        f += (T[i]-mu[i])*(T[i+minLag+5]-mu[i+minLag+5]);
        g += (T[i]-mu[i])*(T[i+minLag+6]-mu[i+minLag+6]);
        h += (T[i]-mu[i])*(T[i+minLag+7]-mu[i+minLag+7]);
    }
*/    
    output[0+j] = a;
    output[1+j] = b;
    output[2+j] = c;
    output[3+j] = d;
    //output[4+j] = e;
    //output[5+j] = f;
    //output[6+j] = g;
    //output[7+j] = h;
    }
}


void singleBlockTest(double* T, double* mu, double* output, int minLag, int n, int m){
    for(int i = 0; i < 8; i++){
        double a  = 0;
        double mL = mu[minLag];
        double nL = mu[i];
        for(int j = 0; j < m; j++){
            a += (T[i+j]-nL)*(T[i+j+minLag]);//-mu[i+j+minLag]);
        }
        output[i] = a;
    }
}

#define statSz 32
void statBufferTest(double* T, double* mu, double* output, int minLag, int n, int m){
    double t[statSz];
    double nL = mu[0];
    double wL = mu[minLag];
    for(int j = 0; j < m; j++){
        for(int i = 0; i < statSz; i++){
            t[i] += (T[j]-nL)*(T[j+minLag+i]-mu[i+j+minLag]);
        }
    }
    for(int i = 0; i < statSz; i++){
        output[i] = t[i];
    }
}


void nonav(double* T, double* mu, double* output, int minLag, int n, int m){
    for(int j = 0; j < 32*m; j+=4){
      /* double a = T[j];
       double b = T[j+1];
       double c = T[j+2];
       double d = T[j+3];
       a += T[j+minLag];
       b += T[j+minLag];
       c += T[j+minLag];
       d += T[j+minLag];
       output[j] = a;
       output[j+1] = b;
       output[j+2] = c;
       output[j+3] = d;*/
       for(int i = 0; i < 4; i++){
           T[i+j] += T[i+j+minLag];
       }
    }
}


#define vsz 4
#ifdef __AVX2__
void avtest(double* T, double* mu, double* output, int minLag, int n, int m){

    double t[statSz];
    for(int j = 0; j < m; j++){
        for(int i = 0; i < statSz; i++){
            __m256d ac2 = _mm256_load_pd(&t[i]);
            __m256d ac0 = _mm256_loadu_pd(&T[j]);
            __m256d ac1 = _mm256_loadu_pd(&T[j+minLag]);
            //__m256d m0  = _mm256_loadu_pd(&mu[j]);
            //__m256d m1  = _mm256_loadu_pd(&mu[j+minLag]);
            //__m256d ac2 = _mm256_sub_pd(ac0,m0);
            //__m256d ac3 = _mm256_sub_pd(ac1,m1);
            ac0 = _mm256_fmadd_pd(ac0,ac1,ac2);
            _mm256_store_pd(&t[i],ac0);
        }
        
    }
    memcpy(output,t,statSz*4);
/*
    for(int i = 0; i < statSz; i+=4){
        
    }
  
 */ 

   // for(int j = 0; j < 32*m; j+=4){
    //    __m256d ac1 = _mm256_loadu_pd(&T[j]);
       /* __m256d ac2 = _mm256_loadu_pd(&T[j+4]);
        __m256d ac3 = _mm256_loadu_pd(&T[j+8]); 
        __m256d ac4 = _mm256_loadu_pd(&T[j+12]);
        __m256d ac5 = _mm256_loadu_pd(&T[j+minLag]);
        __m256d ac6 = _mm256_loadu_pd(&T[j+4+minLag]);
        __m256d ac7 = _mm256_loadu_pd(&T[j+8+minLag]);
        __m256d ac8 = _mm256_loadu_pd(&T[j+12+minLag]);
        */
      //  __m256d ac5 = _mm256_loadu_pd(&T[j+minLag]);
      //  _mm256_storeu_pd(&output[j],_mm256_fmadd_pd(ac1,ac5,ac1));
        /*ac1 = _mm256_fmadd_pd(ac1,ac5,ac1);
        ac2 = _mm256_fmadd_pd(ac2,ac6,ac2);
        ac3 = _mm256_fmadd_pd(ac3,ac7,ac3);
        ac4 = _mm256_fmadd_pd(ac4,ac8,ac4);
        _mm256_storeu_pd(&output[j],ac1);
        _mm256_storeu_pd(&output[j],ac2);       
        _mm256_storeu_pd(&output[j],ac3);
        _mm256_storeu_pd(&output[j],ac4); 
        */
   // }
}

/*
void avtest2(double* T, double* mu, double* output, int minLag, int n, int m){
        __m256d ymm0 = _mm256_broadcast_sd(&T[0]);
        __m256d ymm1 = _mm256_loadu_pd(&T[minLag]);
        __m256d ymm2 = _mm256_loadu_pd(&T[minLag+4]);
        __m256d ymm3 = _mm256_loadu_pd(&T[minLag+8]);
        __m256d ymm4 = _mm256_loadu_pd(&T[minLag+12]);
        __m256d ymm5 = _mm256_loadu_pd(&T[minLag+16]);
        __m256d ymm6 = _mm256_loadu_pd(&T[minLag+20]);
        __m256d ymm7 = _mm256_loadu_pd(&T[minLag+24]);
        __m256d ymm8 = _mm256_loadu_pd(&T[minLag+28]);
        __m256d ymm9 = _mm256_loadu_pd(&T[minLag+32]);
        __m256d ymm10 =_mm256_loadu_pd(&T[minLag+36]);        
        __m256d ymm11 =_mm256_loadu_pd(&T[minLag+40]);
        
    for(int i = 0; i < 1024; i++){ 
        ymm1 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+i]),ymm1);
        ymm2 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+4+i]),ymm2);
        ymm3 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+8+i]),ymm3);
        ymm4 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+12+i]),ymm4);
        ymm5 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+16+i]),ymm5);
        ymm6 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+20+i]),ymm6);
        ymm7 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+24+i]),ymm7);
        ymm8 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+28+i]),ymm8);
        ymm9 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+32+i]),ymm9);
        ymm10= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+36+i]),ymm10);
        ymm11= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&mu[minLag+40+i]),ymm11);
  
    }
    _mm256_storeu_pd(output,ymm1);
    _mm256_storeu_pd(output+4,ymm2);
    _mm256_storeu_pd(output+8,ymm3);
    _mm256_storeu_pd(output+12,ymm4);
    _mm256_storeu_pd(output+16,ymm5);
    _mm256_storeu_pd(output+20,ymm6);
    _mm256_storeu_pd(output+24,ymm7);
    _mm256_storeu_pd(output+28,ymm8);
    _mm256_storeu_pd(output+32,ymm9);
    _mm256_storeu_pd(output+36,ymm10);
    _mm256_storeu_pd(output+40,ymm11);
}*/

void avtest2(double* T, double* mu, double* output, int minLag, int n, int m){
        double t[512];
        __m256d ymm0 = _mm256_broadcast_sd(&T[0]);
        __m256d ymm1 = _mm256_loadu_pd(&T[minLag]);
        __m256d ymm2 = _mm256_loadu_pd(&T[minLag+4]);
        __m256d ymm3 = _mm256_loadu_pd(&T[minLag+8]);
        __m256d ymm4 = _mm256_loadu_pd(&T[minLag+12]);
        __m256d ymm5 = _mm256_loadu_pd(&T[minLag+16]);
     /*   __m256d ymm6 = _mm256_loadu_pd(&T[minLag+20]);
        __m256d ymm7 = _mm256_loadu_pd(&T[minLag+24]);
        __m256d ymm8 = _mm256_loadu_pd(&T[minLag+28]);
        __m256d ymm9 = _mm256_loadu_pd(&T[minLag+32]);
        __m256d ymm10 =_mm256_loadu_pd(&T[minLag+36]);
        __m256d ymm11 =_mm256_loadu_pd(&T[minLag+40]);
*/
        _mm256_storeu_pd(output,ymm1);
        _mm256_storeu_pd(output+4,ymm2);
        _mm256_storeu_pd(output+8,ymm3);
        _mm256_storeu_pd(output+12,ymm4);
        _mm256_storeu_pd(output+16,ymm5);
       /* _mm256_storeu_pd(output+20,ymm6);
        _mm256_storeu_pd(output+24,ymm7);
        _mm256_storeu_pd(output+28,ymm8);
        _mm256_storeu_pd(output+32,ymm9);
        _mm256_storeu_pd(output+36,ymm10);
        _mm256_storeu_pd(output+40,ymm11);
*/    for(int i = 0; i < 4*1024; i++){
        ymm1 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*4]),ymm1);
        ymm2 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*8]),ymm2);
        ymm3 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*12]),ymm3);
        ymm4 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*16]),ymm4);
        ymm5 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*20]),ymm5);
  /*      ymm6 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*24]),ymm6);
        ymm7 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*28]),ymm7);
        ymm8 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*32]),ymm8);
        ymm9 = _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*36]),ymm9);
        ymm10= _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*40]),ymm10);
        ymm11= _mm256_fmadd_pd(ymm0,_mm256_load_pd(&T[i*44]),ymm11);
*/
    
        _mm256_storeu_pd(t,ymm1);
        _mm256_storeu_pd(t+4,ymm2);
        _mm256_storeu_pd(t+8,ymm3);
        _mm256_storeu_pd(t+12,ymm4);
        _mm256_storeu_pd(t+16,ymm5);
 /*       _mm256_storeu_pd(t+20,ymm6);
        _mm256_storeu_pd(t+24,ymm7);
        _mm256_storeu_pd(t+28,ymm8);
        _mm256_storeu_pd(t+32,ymm9);
        _mm256_storeu_pd(t+36,ymm10);
        _mm256_storeu_pd(t+40,ymm11);
  */  }
    memcpy(output,t,512*4);
}



void avtest3(double* T, double* mu, double* sInv, double* output, int minLag, int n, int m){
        //double t[2*512];
     //   __m256d ymm0 = _mm256_broadcast_sd(&T[0]);
        __m256d ymm1 = _mm256_loadu_pd(&T[minLag]);
        __m256d ymm2 = _mm256_loadu_pd(&T[minLag+4]);
        __m256d ymm3 = _mm256_loadu_pd(&T[minLag+8]);
        __m256d ymm4 = _mm256_loadu_pd(&T[minLag+12]);
        __m256d ymm5 = _mm256_loadu_pd(&T[minLag+16]);
        __m256d ymm6 = _mm256_loadu_pd(&T[minLag+20]);
        __m256d ymm7 = _mm256_loadu_pd(&T[minLag+24]);
        __m256d ymm8 = _mm256_loadu_pd(&T[minLag+28]);
        __m256d ymm9 = _mm256_loadu_pd(&T[minLag+32]);
        __m256d ymm10 =_mm256_loadu_pd(&T[minLag+36]);
        __m256d ymm11 =_mm256_loadu_pd(&T[minLag+40]);
        __m256d ymm12 =_mm256_loadu_pd(&T[minLag+48]);
        __m256d ymm13 =_mm256_loadu_pd(&T[minLag+52]); 
        __m256d ymm14 =_mm256_loadu_pd(&T[minLag+56]);

    for(int j = 0; j < 16; j++){
    for(int i = 0; i < 1024; i++){
        
        __m256d ymm0 = _mm256_broadcast_sd(&T[i]);
 
        ymm1 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*4]),ymm1);
        ymm2 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*8]),ymm2);
        ymm3 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*12]),ymm3);
        ymm4 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*16]),ymm4);
        ymm5 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*20]),ymm5);
        ymm6 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*24]),ymm6);
        ymm7 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*28]),ymm7);
        ymm8 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*32]),ymm8);
        ymm9 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*36]),ymm9);
        ymm10= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*40]),ymm10);
        ymm11= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*44]),ymm11);
        ymm12= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*48]),ymm12);
        ymm13= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*52]),ymm13);
        ymm14= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*56]),ymm14);

        ymm1 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*4]),ymm1);
        ymm2 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*8]),ymm2);
        ymm3 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*12]),ymm3);
        ymm4 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*16]),ymm4);
        ymm5 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*20]),ymm5);
        ymm6 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*24]),ymm6);
        ymm7 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*28]),ymm7);
        ymm8 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*32]),ymm8);
        ymm9 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*36]),ymm9);
        ymm10= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*40]),ymm10);
        ymm11= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*44]),ymm11);
        ymm12= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*48]),ymm12);
        ymm13= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*52]),ymm13);
        ymm14= _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*56]),ymm14);

        
        ymm1 = _mm256_max_pd(ymm1,ymm2);
        ymm3 = _mm256_max_pd(ymm3,ymm4);
        ymm5 = _mm256_max_pd(ymm5,ymm6);
        ymm7 = _mm256_max_pd(ymm7,ymm8);
        ymm9 = _mm256_max_pd(ymm9,ymm10);
        ymm11 = _mm256_max_pd(ymm11,ymm12);
        ymm13 = _mm256_max_pd(ymm13,ymm14);

        ymm2 = _mm256_max_pd(ymm1,ymm3);
        ymm4 = _mm256_max_pd(ymm5,ymm7);
        ymm6 = _mm256_max_pd(ymm9,ymm11);
        ymm1 = _mm256_max_pd(ymm2,ymm4);
        ymm3 = _mm256_max_pd(ymm6,ymm13);
        

        _mm256_storeu_pd(output+4*j+4*i,_mm256_max_pd(ymm1,ymm3));
         
        /*_mm256_storeu_pd(t+j*4,ymm1);
        _mm256_storeu_pd(t+j*8,ymm2);
        _mm256_storeu_pd(t+j*12,ymm3);
        _mm256_storeu_pd(t+j*16,ymm4);
        _mm256_storeu_pd(t+j*20,ymm5);
        _mm256_storeu_pd(t+j*24,ymm6);
        _mm256_storeu_pd(t+j*28,ymm7);
        _mm256_storeu_pd(t+j*32,ymm8);
        _mm256_storeu_pd(t+j*36,ymm9);
        _mm256_storeu_pd(t+j*40,ymm10);
        _mm256_storeu_pd(t+j*44,ymm11);
        */
    }
    }
    //memcpy(output,t,512*4);
}


void avtest4(double* T, double* mu, double* output, int minLag, int n, int m){
        //double t[2*512];
     //   __m256d ymm0 = _mm256_broadcast_sd(&T[0]);
        __m256d ymm1 = _mm256_loadu_pd(&T[minLag]);
        __m256d ymm2 = _mm256_loadu_pd(&T[minLag+4]);
        __m256d ymm3 = _mm256_loadu_pd(&T[minLag+8]);
        __m256d ymm4 = _mm256_loadu_pd(&T[minLag+12]);
        __m256d ymm5 = _mm256_loadu_pd(&T[minLag+16]);
        __m256d ymm6 = _mm256_loadu_pd(&T[minLag+20]);
        __m256d ymm7 = _mm256_loadu_pd(&T[minLag+24]);
        __m256d ymm8 = _mm256_loadu_pd(&T[minLag+28]);
        __m256d ymm9 = _mm256_loadu_pd(&T[minLag+32]);
        __m256d ymm10 =_mm256_loadu_pd(&T[minLag+36]);
        __m256d ymm11 =_mm256_loadu_pd(&T[minLag+40]);
        __m256d ymm12 =_mm256_loadu_pd(&T[minLag+48]);
        __m256d ymm13 =_mm256_loadu_pd(&T[minLag+52]); 
        __m256d ymm14 =_mm256_loadu_pd(&T[minLag+56]);

    for(int j = 0; j < 16; j++){
    for(int i = 0; i < 1024; i++){
        
        __m256d ymm0 = _mm256_broadcast_sd(&T[i]);
 
        _mm256_storeu_pd(output+4*j+4*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*4]),ymm1));
        _mm256_storeu_pd(output+8*j+8*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*8]),ymm2));
        _mm256_storeu_pd(output+12*j+12*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*12]),ymm3));
        _mm256_storeu_pd(output+16*j+16*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*16]),ymm4));
        _mm256_storeu_pd(output+20*j+20*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*20]),ymm5));
        _mm256_storeu_pd(output+24*j+24*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*24]),ymm6));
        _mm256_storeu_pd(output+28*j+28*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*28]),ymm7));
        _mm256_storeu_pd(output+32*j+32*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*32]),ymm8));
        _mm256_storeu_pd(output+36*j+36*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*36]),ymm9));
        _mm256_storeu_pd(output+40*j+40*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*40]),ymm10));
        _mm256_storeu_pd(output+44*j+44*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*44]),ymm11));
        _mm256_storeu_pd(output+48*j+48*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*48]),ymm12));
        _mm256_storeu_pd(output+52*j+52*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*52]),ymm13));
        _mm256_storeu_pd(output+56*j+56*i,_mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[i*56]),ymm14));

    }
    }
    //memcpy(output,t,512*4);
}

void avtest5(const double* T,const  double* mu, const double* s, double* output, int minLag, int n, int m){
    for(int i = 0; i < 1024; i++){
        __m256d ymm0 = _mm256_broadcast_sd(&T[i]);
        __m256d ymm4 = _mm256_loadu_pd(&T[minLag+4*i]);
        __m256d ymm2 = _mm256_loadu_pd(&s[4*i]);
        __m256d ymm1 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[minLag+4*i]),ymm4);
        __m256d ymm5 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[minLag+4*i+4]),ymm4);
        ymm1 = _mm256_mul_pd(ymm1,ymm2);
        for(int j = 1; j < 16; j++){
            for(int k = j; k < j+14; k++){
               __m256d ymm4 = _mm256_loadu_pd(&T[minLag+4*(i+k)]);
               __m256d ymm2 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[minLag+4*i+4*k]),ymm4);
               ymm5 = _mm256_fmadd_pd(ymm0,_mm256_loadu_pd(&T[minLag+8*i+4*k]),ymm2);
               __m256d ymm3 = _mm256_loadu_pd(&s[4*(i+k)]);
               ymm2 = _mm256_mul_pd(ymm5,ymm3);
               ymm1 = _mm256_max_pd(ymm1,ymm2);
            }
        }
        _mm256_storeu_pd(output+4*i,ymm1);
    }
}

#endif




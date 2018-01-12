#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<math.h>
#include<omp.h>

//extern void  accumtest4_7_9_loop(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n, int m);
extern void  accumtest4_7_9(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n, int m);
extern void simd_initstats(double* a, double* mu, double* s, int n, int w);
//extern void winmeanalt(double* a, double* mu, double* s, int n, int w);
//extern void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n,int m);

void writeDoubles(const char* name, double* t, int n){
   FILE* f = fopen(name,"w");
   if(f == NULL){
      perror("fopen");  
      exit(1);
   }
   for(int i = 0; i < n; i++){
      fprintf(f,"%.15lf\n",t[i]);
   }
   fclose(f);
}

void writeInts(const char* name, int* t, int n){
   FILE* f = fopen(name,"w");
   if(f == NULL){
      perror("fopen");  
      exit(1);
   }
   for(int i = 0; i < n; i++){
      fprintf(f,"%d\n",t[i]);
   }
   fclose(f);
}


int main(int argc, char* argv[]){
    if(argc < 4){
        printf("check input arguments\n");
        exit(0);
    }
    int n = atoi(argv[2]);
    FILE* f = fopen(argv[1],"r");
    int m = atoi(argv[3]);
    
    if(f == NULL){
        perror("fopen");
        exit(1);
    }
    double* T = NULL;
    posix_memalign((void**)&T,64,n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }    
    printf("check1\n");
    fclose(f);

    double* mu = NULL;
    double* buffer = NULL;
    double* sigmaInv = NULL;
    double* dX = NULL;
    double* dF = NULL;
    double* bufA = NULL;
    double* bufB = NULL;
    long*   bufferI = NULL;
    posix_memalign((void**)&mu,64,8*n*sizeof(double));
    posix_memalign((void**)&sigmaInv,64,8*n*sizeof(double));
/*    posix_memalign((void**)&buffer,64,8*n*sizeof(double));
    posix_memalign((void**)&dX,64,8*n*sizeof(double));
    posix_memalign((void**)&dF,64,8*n*sizeof(double));
    posix_memalign((void**)&bufferI,64,8*n*sizeof(long));
*/
 /*   for(int i = 0; i < n; i++){
        buffer[i] = -1.0;
    }
    clock_t t1 = clock();
    clock_t t2 = clock();
*/  printf("%lf %lf %lf %lf\n",T[n-1],dX[n-m-1],dF[n-m-1],sigmaInv[n-m-1]);
    printf("done\n");
    printf("test:  %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    writeDoubles("meancheck.csv",mu,n-m+1);
    writeDoubles("scheck.csv", sigmaInv,n-m+1); 
    free(T);
//    free(buffer);
//    free(bufferI);
    free(mu);
    free(sigmaInv);
    //free(dX);
    //free(dF);
    return 0;
}

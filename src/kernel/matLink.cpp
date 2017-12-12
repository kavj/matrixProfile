#define _POSIX_C_SOURCE 200809L
#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include<math.h>
#include<omp.h>



/* Use wc -l <filename> to get the number of lines in a single column csv. Pass it as argument 2 here.*/

/*extern tsdesc* sc_init(double* T, int n, int m);
extern void sc_destroy(tsdesc* T);
extern matrixProfileObj* mp_init(int n, int m);
extern void mp_destroy(matrixProfileObj* mp); 
*/

extern void accumTest4(double* Cxy, double* dX, double* Sxy, double* output, double* fw, double* bw, int n, int m);

void  accumTest4_10(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n);
extern void  accumTest4_7_16(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);
extern void  accumTest4_7_12(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);
extern void  accumTest4_7_18(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);
extern void  accumTest4_7_16alt(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);


extern void  accumTest4_7_11_alt(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);
extern void  accumTest4_7_11(double* cxy, double* dx, double* df, double*s, double* output, long* outputi,int n);

extern void  accumTest4_7_2(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n);
extern void  accumTest4_7_9(double* Cxy, double* dX, double* dF, double*s, double* output, long* outputI,int n);
extern void  shortconv(double* a, double* b, double* c, int n);


static inline double shiftMean(double mu, double a, double b, double m){
    return mu + (a-b)/m;
}

static inline double shiftSSS(double s, double mu, double muprev, double x, double xprev){
    return s + ((x - mu) * (x - muprev) - (xprev - mu) * (xprev - muprev));
}

static inline double shiftSXY(double x, double y, double xprev, double yprev, double mux, double muyprev){
    return (x - mux)*(y - muyprev) - (xprev - mux)*(yprev - muyprev);
}

static inline double initSS(const double* T, double M, int m, int offset){
    double s = 0;
    for(int i = offset; i < offset + m; i++){
        s += (T[i] - M) * (T[i] - M);
    }
    return s;
}


static inline double initMean(const double* T, int m, int offset){
    double M = T[offset];
    for(int i = offset+1; i < offset+m; i++){
        M += (T[i] - M)/(i-offset+1);
    }
    return M;
}

/* This uses a variation on Welford's method for sample variance to compute the mean and standard deviation. 
 *  * It resets the summation once the term under consideration does not share any terms with the last exact summation. */
/* This should check for exceptions*/
/* I should make this clear somewhere that this isn't a full standard deviation function. It could be refactored into one,
 *  * but it should preserve the m factor difference given that this is used to cancel another similar factor */
void winmeansig(const double* T, double* mu, double* sigma, int n, int m,int normConstant){
    int alignedBound = (n-m+1) - (n-m+1) % m;
    for(int i = 0; i < alignedBound; i += m){
        double M    = initMean(T,m,i);
        double s    = initSS(T,M,m,i);
        mu[i]       = M;
        sigma[i]    = sqrt(s/normConstant);

        for(int j = i+1; j < i+m; j++){
            double Mprev = M;
            M        = shiftMean(M,T[j+m-1],T[j-1],m);
            s        = shiftSSS(s,M,Mprev,T[j+m-1],T[j-1]);
            mu[j]    = M;
            sigma[j] = sqrt(s/normConstant);
        }
    }

    /* compute unaligned portion */
    double M               = initMean(T,m,alignedBound);
    double s               = initSS(T,M,m,alignedBound);
    mu[alignedBound]       = M;
    sigma[alignedBound]    = s;
    for(int i = alignedBound; i < n-m+1; i++){
        sigma[i] = sqrt(s/normConstant);
    }
}


void initDXDF(double* x,double* mu, double* dF, double* dX,int n, int m){
    for(int i = 0; i < n-m; i++){
        dX[i] = x[i+m]-x[i];
        dF[i] = (x[i+m]-mu[i+1]) + (x[i]-mu[i]);
    }
}


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
    //double* T = malloc(n*sizeof(double));
    double* T = NULL;
    posix_memalign((void**)&T,64,n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }    
    printf("check1\n");
    fclose(f);
    writeDoubles("sanitycheckT.csv",T,n);

    double* mu = NULL;
    double* buffer = NULL;
//    double* bufferB = NULL;
    double* sigmaInv = NULL;
    double* dX = NULL;
    double* dF = NULL;
    long*   bufferI = NULL;
    posix_memalign((void**)&mu,64,8*n*sizeof(double));
    posix_memalign((void**)&buffer,64,8*n*sizeof(double));
//    posix_memalign((void**)&bufferB,64,4*n*sizeof(double));
    posix_memalign((void**)&sigmaInv,64,8*n*sizeof(double));
    posix_memalign((void**)&dX,64,8*n*sizeof(double));
    posix_memalign((void**)&dF,64,8*n*sizeof(double));
    posix_memalign((void**)&bufferI,64,16*n*sizeof(long));
    for(int i = 0; i < n; i++){
        buffer[i] = -2*m;
    }
    clock_t t1 = clock();
    winmeansig(T,mu,sigmaInv,n,m,m);
    initDXDF(T, mu, dF, dX, n, m);
    clock_t t2 = clock();
    printf("%lf %lf %lf %lf\n",T[n-1],dX[n-m-1],dF[n-m-1],sigmaInv[n-m-1]);
    clock_t t3 = clock();//omp_get_wtime();
    accumTest4_7_11(T,dX,dF,sigmaInv,buffer,bufferI,n);
    clock_t t4 = clock();//omp_get_wtime();
    printf("done\n");
    printf("test:  %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    printf("test:  %lf\n",(double)(t4-t3)/CLOCKS_PER_SEC);
    free(T);
    free(buffer);
 //   free(bufferB);
    free(bufferI);
    free(mu);
    free(sigmaInv);
    free(dX);
    free(dF);
    return 0;
}

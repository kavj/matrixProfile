#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>
#include "scr.h"


/* Use wc -l <filename> to get the number of lines in a single column csv. Pass it as argument 2 here.*/

/*extern tsdesc* sc_init(double* T, int n, int m);
extern void sc_destroy(tsdesc* T);
extern matrixProfileObj* mp_init(int n, int m);
extern void mp_destroy(matrixProfileObj* mp); 
*/

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
    if(argc < 3){
        printf("check input arguments\n");
        exit(0);
    }
    int n = atoi(argv[2]);
    FILE* f = fopen(argv[1],"r");
    int m = 250;
    
    if(f == NULL){
        perror("fopen");
        exit(1);
    }
    double* T = malloc(n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }    
    fclose(f);
    writeDoubles("sanitycheckT.csv",T,n);
    tsdesc* t = sc_init(T,n,m);
    matrixProfileObj* matp = mp_init(n,m);
    clock_t t1 = clock();
    winmeansig(t->T,t->mu,t->sigmaInv,n,m);
    scBlockSolver(t,matp);
    //corrToDist(matp->mp,n,m);
    clock_t t2 = clock();
    printf("time: %lf\n",(double)(t2-t1)/CLOCKS_PER_SEC);
    writeDoubles("penguin_mp.csv",matp->mp,n-m+1);
    writeInts("penguin_mpI.csv",matp->mpI,n-m+1);
    writeDoubles("penguin_mu.csv", t->mu,n-m+1);
    writeDoubles("penguin_sigma.csv",t->sigmaInv,n-m+1);
    sc_destroy(t);
    mp_destroy(matp);
    free(T);
    return 0;
}

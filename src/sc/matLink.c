#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include "scr.h"

/* Use wc -l <filename> to get the number of lines in a single column csv. Pass it as argument 2 here.*/

/*extern tsdesc* sc_init(double* T, int n, int m);
extern void sc_destroy(tsdesc* T);
extern matrixProfileObj* mp_init(int n, int m);
extern void mp_destroy(matrixProfileObj* mp); 
*/

void writeFile(const char* name, double* t, int n){
    FILE* f = fopen(name,"w");
    if(f == NULL){
        perror("fopen");  
        exit(1);
    }
    for(int i = 0; i < n; i++){
        fprintf(f,"%lf\n",t[i]);
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
    }
    double* T = malloc(n*sizeof(double));
    for(int i = 0; i < n; i++){
        fscanf(f,"%lf\n",&T[i]);
    }    
    fclose(f);
    tsdesc* t = sc_init(T,n,m);
    matrixProfileObj* matp = mp_init(n,m);
    winmeansig(t->T,t->mu,t->sigmaInv,n,m);
    scBlockSolver(t,matp);

    writeFile("penguin_mp.csv",matp->mp);
    writeFile("penguin_mpI.csv",matp->mpI);
    writeFile("penguin_mu.csv", t->mu);
    writeFile("penguin_sigma.csv",t->sigmaInv);
    sc_destroy(t);
    free(T);
    return 0;
}

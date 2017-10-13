#include<stdio.h>
#include<stdlib.h>



static void compMeanSS(const double* T, double* mu, double* s, int n, int m){
    for(int i = 0; i < n-m; i++){
        mu[i] = T[i];
        for(int j = 0; j < m; j++){
            mu[i] += (T[i+j]-mu[i])/j;
        }
        s[i] = 0;
        for(int j = 0; j < m; j++){
            s[i] += (T[i+j] - mu[i])*(T[i+j]-mu[i]);
        }
    }
}

static void twoPass(const double* T, const double* mu, const double* s, double* c, int n , int m, int lag){
    for(int i = 0; i < m; i++){
        double d = 0;
        double mup = mu[i];
        double mulag = mu[i+lag];
        for(int j = 0; j < m; j++){
            d += (T[i+j] - mup)*(T[i+j+lag] - mulag);
        }
        c[i] = d*s[i]*s[i+lag];
    }
}
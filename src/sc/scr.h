#ifndef SCR
#define SCR

struct desc{
    double* T;
    double* mu;
    double* sigmaInv;
    int n;
    int m;
    int chsz;
};

struct mpObj{
    double* mp;
    int* mpI;
};


typedef struct desc tsdesc;
typedef struct mpObj matrixProfileObj;


int iterCnt(int* x, int xOffs, int n, int m, double perc);

void winmeansig(const double* T, double* mu, double* sigma, int n, int m);

void mpSelf(const double* T, const double* mu, const double* sigma, double* mp, int* mpI, int n, int m, int lag); 

tsdesc* sc_init(double* T, int n, int m);

matrixProfileObj* mp_init(int n, int m);

void sc_destroy(tsdesc* t);

void mp_destroy(matrixProfileObj* matp);



#endif

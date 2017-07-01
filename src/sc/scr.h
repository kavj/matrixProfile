#ifndef SCR
#define SCR

truct mpStats{
    const double* T;
    double* mu;
    double* sigma;
    int n;
    int m;
};

struct mpObj{
    double* mp;
    int* mpI;
};


typedef struct mpStats tsdesc;
typedef struct mpObj matrixProfileObj;


int iterCnt(int* x, int xOffs, int n, double perc);

void winmeansig(const double* T, double* mu, double* sigma, int n, int m);

void mpSelf(const double* T, const double* mu, const double* sigma, double* mp, 

tsdesc* sc_init(const double* T, int n, int m);

matrixProfileObj* mp_init(int n, int m);

void sc_destroy(tsdesc* t);

void mp_destroy(matrixProfileObj* matp);



#endif

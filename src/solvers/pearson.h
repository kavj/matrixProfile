
typedef struct profile{
   double* matrixProfile;
   long long* matrixProfileIndex;
   long long len;
} mprofile;


mprofile pearson_pauto_reduc(const double* ts, long long len, long long minlag, long long sublen); 



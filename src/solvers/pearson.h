
typedef struct profile{
   double* matrixProfile;
   long long* matrixProfileIndex;
} mprofile;


mprofile pearson_pauto_reduc(const double* ts, long long len, long long minlag, long long sublen); 


//Todo: Figure out what to do with error handling
/*namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};*/

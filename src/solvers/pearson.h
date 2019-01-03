


int nautocorr_reduc(double* ts, double*  mp, long long* mpi, long long minlag, long long sublen);

void pearson2zned(double* mp, long long len, long long sublen);


namespace errs{
   const int bad_input = -1;
   const int mem_error = -1;
   const int none = 0;

};

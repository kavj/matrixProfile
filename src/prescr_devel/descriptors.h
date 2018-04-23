


// Todo: figure out a reasonable way to turn factoor out separate back ends for mkl vs quasi Cook-Toom filtration


struct qbuf{
   double* qbuf;  // buffer for multiple queries
   double* qcov;  // covariance of optimal query match, multiple copies are used in the case of multiple queries
   double* qcorr; // need to know this
//   double* qmu;   // may be helpful when 
   int* qbase;    // stores the base index of each query with respect to the time series from which it was sampled
   int* qmatch;   // nearest neighbor index for each query
   int qtotal;
   int qbufcount;
   int qlen;
   int qstride;
   int qbufstride;
};


struct corr_desc{
   double* xcorr;
   int* xcorrind;
   int ccount;     //number of elements in correlation vector, in our cases (we truncate aliased portions) length(time_series) - length(query) + 1
   int clen;       //length of reduced cross correlation vector
   int cseclen;    //length of each correlation block section
   int cbufstride; //stride with respect to correlation buffer, a stride in memory accounts for fringe portions and elements used to pad memory alignment
   
};





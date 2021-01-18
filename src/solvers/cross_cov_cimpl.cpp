#include "mex.h"
#include "matrix.h"
#include "cross_cov.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[] ){
    // cv, r_bwd, c_bwd, r_fwd, c_fwd, invn_r, invn_c

    if(nrhs != 7){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "(required signature is (ts, mu, cmpseq).");
    }
    else if (nlhs != 1){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumOutputs", "This is meant to output cv.");
    }
    
    // This shouldn't receive an entire time series on self comparisons.
    // cmpseq is whatever sequence is uniform across all initial comparisons
    // ts contains all sequences that are allowed to compare to cmpseq. For the "self-join"
    // logic, this means ts(minsep:end), mu(minsep:end)
    double* ts = mxGetDoubles(prhs[0]);
    double* mu = mxGetDoubles(prhs[1]);
    double* cmpseq = mxGetDoubles(prhs[2]);
    long long tslen = mxGetNumberOfElements(prhs[0]);
    long long dcount = mxGetNumberOfElements(prhs[1]);
    long long sublen = mxGetNumberOfElements(prhs[2]);
    if(tslen - sublen + 1 != dcount){
        mexErrMsgIdAndTxt( "MATLAB:invalidlength", "Check inputs. Length of mu does not match subsequence count in ts.");
    }

    mxArray* cv_ = mxCreateDoubleMatrix(dcount, 1, mxREAL);
    double* cv = mxGetDoubles(cv_);

    // co-moments contain comparison between index 0 and index minsep....subcount
 
    crosscov(cv, ts, mu, cmpseq, dcount, sublen);

    plhs[0] = cv_;
    return;
}


#include "mex.h"
#include "matrix.h"
#include "pearson_split.h"

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray* prhs[] ){
    // cv, r_bwd, c_bwd, r_fwd, c_fwd, invn_r, invn_c

    if(nrhs != 7){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumInputs", "(required signature is (cv, dr_bwd, dc_bwd, dr_fwd, dc_fwd, invn_r, invn_c, mp_r, mp_c).");
    }
    else if (nlhs != 2){
        mexErrMsgIdAndTxt( "MATLAB:invalidNumOutputs", "This is meant to output mp_rows, mp_cols.");
    }
    
    // look up check for scalar int on third 
    
    double* cv = mxGetDoubles(prhs[0]);
    double* dr_bwd = mxGetDoubles(prhs[1]);
    double* dc_bwd = mxGetDoubles(prhs[2]);
    double* dr_fwd = mxGetDoubles(prhs[3]);
    double* dc_fwd = mxGetDoubles(prhs[4]);
    double* invn_r = mxGetDoubles(prhs[5]);
    double* invn_c = mxGetDoubles(prhs[6]);
    long long dcount = mxGetNumberOfElements(prhs[0]);

    mxArray* mp_r_ = mxCreateDoubleMatrix(dcount, 1, mxREAL);
    mxArray* mp_c_ = mxCreateDoubleMatrix(dcount, 1, mxREAL);

    double* mp_r = mxGetDoubles(mp_r_);
    double* mp_c = mxGetDoubles(mp_c_);

    // co-moments contain comparison between index 0 and index minsep....subcount
 
    compute_self_mp(cv, mp_r, mp_c, dr_bwd, dc_bwd, dr_fwd, dc_fwd, invn_r, invn_c, dcount);

    plhs[0] = mp_r_;
    plhs[1] = mp_c_;
    return;
}


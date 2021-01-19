#include<algorithm>
#include "mex.h"
#include "matrix.h"
#include "pearson_split.h"


extern void mpx_ref(double* mpr, double* mpc, double* cv, double* dr_bwd, double* dc_bwd, double* dr_fwd, double* dc_fwd, double* invnr, double* invnc, int dcount);

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

    std::fill(mp_r, mp_r+dcount, -1.0);
    std::fill(mp_c, mp_c+dcount, -1.0);

    // co-moments contain comparison between index 0 and index minsep....subcount

    long long rbct = mxGetNumberOfElements(prhs[1]);
    long long cbct = mxGetNumberOfElements(prhs[2]);
    long long rfct = mxGetNumberOfElements(prhs[3]);
    long long cfct = mxGetNumberOfElements(prhs[4]);
    long long invnrct = mxGetNumberOfElements(prhs[5]);
    long long invncct = mxGetNumberOfElements(prhs[6]);
    long long mprct = mxGetNumberOfElements(mp_r_);
    long long mpcct = mxGetNumberOfElements(mp_c_);

    mpx_ref(mp_r, mp_c, cv, dr_bwd, dc_bwd, dr_fwd, dc_fwd, invn_r, invn_c, dcount);
    mexPrintf("dcount:%d rbwd:%d cbwd:%d rfwd:%d cfwd:%d invnrct:%d invncct:%d mprct:%d mpcct%d\n", dcount, rbct, cbct, rfct, cfct, invnrct, invncct, mprct, mpcct);
    
    //compute_self_mp(cv, mp_r, mp_c, dr_bwd, dc_bwd, dr_fwd, dc_fwd, invn_r, invn_c, dcount);

    plhs[0] = mp_r_;
    plhs[1] = mp_c_;
    return;
}


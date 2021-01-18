#ifndef PEARSON_SPLIT
#define PEARSON_SPLIT

void compute_self_mp(double* __restrict cv,
                double* __restrict mpr,
                double* __restrict mpc,
                double* __restrict dr_bwd,
                double* __restrict dc_bwd,
                double* __restrict dr_fwd,
                double* __restrict dc_fwd,
                double* __restrict invnr,
                double* __restrict invnc,
                int dcount);
#endif

#include "mex.h"
#include <iostream>
#include <vector>
#include <cmath>

void
mexFunction (int nlhs, mxArray *plhs[],
             int nrhs, const mxArray *prhs[])  //rhs: J, beta_max, interaction_mat, T, SEED
{
    double J = *mxGetPr(prhs[0]);
    double beta_max = *mxGetPr(prhs[1]);
    bool *interaction_mat = (bool *)mxGetPr(prhs[2]);
    int T = (int)(*(mxGetPr(prhs[3])));
    mwSize n_cells = mxGetM(prhs[2]);
    mwSize dims[2] = {n_cells, (mwSize)1};
    //srand( (unsigned) (*(mxGetPr(prhs[4])));
    srand((unsigned)*mxGetPr(prhs[4]));
    plhs[0] = mxCreateNumericArray(2, dims, mxGetClassID(prhs[2]), mxREAL);
    bool *q = (bool *)mxGetPr(plhs[0]);
    for(int i=0; i<n_cells; i++) q[i] = (rand() % 2 == 1);
    
    //std::cout << "starting" << std::endl;
    
    mxArray *i_times_q_arr = mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxREAL);
    double *i_times_q = mxGetPr(i_times_q_arr);
    for(int i=0; i<n_cells; i++) {
        i_times_q[i] = 0;
        for(int j=0; j<n_cells; j++) {
            i_times_q[i] += (q[j] * interaction_mat[j + n_cells*i]);
        }
    }
    //std::cout << "starting t loop, J=" << J << " beta_max=" << beta_max<< " T=" << T << std::endl;
    for(int t = 1; t <= T; t++) {
        double beta = ((double)t)/((double)T) * beta_max;
        //std::cout << "beta " << beta << std::endl;
        for (int i=0; i<n_cells; i++) {
            int delta_q = 1 - 2*q[i];
            double delta_H = delta_q*(J*i_times_q[i]-1);//-delta_q + J*delta_q*i_times_q[i];
            //std::cout << "i=" << i << " delta_q=" << delta_q << " delta_H=" << delta_H << std::endl;
            if((delta_H < 0) || ((((double)rand()) / ((double)(RAND_MAX))) < exp(-beta*delta_H))) {
                q[i] = 1-q[i];
                for(int j=0; j<n_cells; j++) {
                    i_times_q[j] += delta_q*interaction_mat[j + n_cells*i];
                }
            }
        }
    }
}
/*
    //double *r = mxGetPr(prhs[0]);
    //double epsilon = *mxGetPr(prhs[1]);
    //double sigma = *mxGetPr(prhs[2]);
    //double m = *mxGetPr(prhs[3]);
    //double L = *mxGetPr(prhs[4]);
    //int dims[2] = {mxGetM(prhs[0]), mxGetN(prhs[0])};
    //int N = dims[0];
    //plhs[0] = mxCreateNumericArray(2, dims, mxGetClassID(prhs[0]), mxREAL);
    //double *a = mxGetPr(plhs[0]);
    for(int i=0; i<N; i++){ a[i] = 0; a[N+i] = 0;}

    double sigma6 = (sigma * sigma * sigma * sigma * sigma * sigma), sigma12 = sigma6*sigma6;

    for(int i=0; i<N; i++){
        double xi = fmod(r[i],L), yi = fmod(r[N+i],L);
        for(int j=0; j<N; j++){
            if(i==j) continue;
            //for a given pair, i and j
            double xj0 = fmod(r[j],L), yj0 = fmod(r[N+j],L);
            double xj[9], yj[9];
            xj[0] = xj0; xj[1] = xj0; xj[2] = xj0; xj[3] = xj0+L; xj[4] = xj0-L; xj[5] = xj[4]; xj[6] = xj[3]; xj[7] = xj[4]; xj[8] = xj[3];
            yj[0] = yj0; yj[1] = yj0 + L; yj[2] = yj0 - L; yj[3] = yj0; yj[4] = yj0; yj[5] = yj[1]; yj[6] = yj[1]; yj[7] = yj[2]; yj[8] = yj[2];

            for(int k = 0; k < 9; k++){
                double diffx = xi-xj[k], diffy = yi-yj[k];
                double r2 = diffx*diffx + diffy*diffy;
                double tempf = epsilon * (12 * sigma12/(r2*r2*r2) - 6 * sigma6)/(r2*r2*r2*r2);
                //force on i
                a[i] += tempf * diffx/m; a[N+i] += tempf * diffy/m;
                //force on j
                //a[j] += -a[i]; a[N+j] += -a[N+i];
            }
        }
    }
}
*/
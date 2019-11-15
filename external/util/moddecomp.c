#include <stdio.h>
#include <stdlib.h>
#include "mex.h"
#include "matrix.h"

void decomp(int q, int *pM, int *pD, int n){
  int d, r;
  if (q > 0){
    d = pM[n];
    r = q % d;
    q = (q-r)/d;
    pD[n] = r;
    decomp(q,pM,pD,++n);
  }
}
#ifdef STANDALONE
int main(int argc, char **argv){
  int i, N;
  int *pM, *pD, nM;
  
  N = (int)atof(argv[1]);

  nM = argc-2;

  pM = calloc(nM, sizeof(int));
  pD = calloc(nM, sizeof(int));

  for (i = 0; i<nM; i++)
    pM[i] = (int)atof(argv[i+2]);

  decomp(N, pM, pD, 0);

  free(pM);
  free(pD);

  return 0;
}

#else

void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs [] ){
  int i, n = 0;
  int N, *pM, *pD; 
  double *pfM, *pfD;
  int nM, nD;

  if (nrhs != 2){
    mexErrMsgTxt("Usage: moddecomp n M.\n");
  }

  N  = (int)*(mxGetPr(prhs[0]));  
  nM = mxGetN(prhs[1]);
  pfM = mxGetPr(prhs[1]);
  
  pM = mxCalloc(nM, sizeof(int));
  for (i=0; i<nM; i++)
    pM[i] = (int)pfM[i];

  pD = mxCalloc(nM, sizeof(int));

  decomp(N, pM, pD, 0);

  plhs[0] = mxCreateDoubleMatrix(1, nM, mxREAL);
  pfD = mxGetPr(plhs[0]);
  for (i = 0; i<nM; i++)
    pfD[i] = (double)pD[i];

  mxFree(pM);
  mxFree(pD);

}
  
  
#endif

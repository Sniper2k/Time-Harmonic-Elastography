#include "mex.h"

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double lambda;             /* input scalar */
    double q_abs;
    double sqQ;
    double *q;                /* 1xN input matrix */
    size_t nx;             /* size of matrix */
    size_t nxx;
    size_t ny;             /* size of matrix */
    double *p;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=1) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }    
    
    q = mxGetPr(prhs[0]);
    lambda = mxGetScalar(prhs[1]);

    /* get dimensions of the input matrix */
    nx = mxGetM(prhs[0]);
    ny = mxGetN(prhs[0]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nx,(mwSize)ny,mxREAL);    
    plhs[1] = mxCreateDoubleMatrix((mwSize)nxx,(mwSize)ny,mxREAL);

    /* get a pointer to the real data in the output matrix */
    p = mxGetPr(plhs[0]);
    nxx = nx/4;
    
    mwSize i;
    mwSize j;
    for (i=0; i<nxx; i++) {
        for (j=0; j<ny; j++) {
            q_abs = q[i+j*nx]*q[i+j*nx] + q[i+nxx+j*nx]*q[i+nxx+j*nx] + q[i+2*nxx+j*nx]*q[i+2*nxx+j*nx] + q[i+3*nxx+j*nx]*q[i+3*nxx+j*nx];
            if (1 < lambda*lambda / q_abs) {
                p[i+j*nx]       = q[i+j*nx];
                p[i+nxx+j*nx]   = q[i+nxx+j*nx];
                p[i+2*nxx+j*nx] = q[i+2*nxx+j*nx];
                p[i+3*nxx+j*nx] = q[i+3*nxx+j*nx];
            } else {
                sqQ = lambda / sqrt(q_abs);
                p[i+j*nx]       = q[i+j*nx] * sqQ;
                p[i+nxx+j*nx]   = q[i+nxx+j*nx] * sqQ;
                p[i+2*nxx+j*nx] = q[i+2*nxx+j*nx] * sqQ;
                p[i+3*nxx+j*nx] = q[i+3*nxx+j*nx] * sqQ;
            }
        }
    }    
}
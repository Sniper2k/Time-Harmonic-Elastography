#include "mex.h"

double SoftShrink (double a, double param) {
    if (a > param) {
        return a - param;
    } else if (a < -param) {
        return a + param;
    } else {
        return 0;
    }
}

/* The gateway function */
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[])
{
    double alpha;             /* input scalar */
    double *A;                /* 1xN input matrix */
    double *B;                /* 1xN input matrix */
    double *c;                /* 1xN input matrix */
    double *d_1;              /* 1xN input matrix */
    double *d_2;              /* 1xN input matrix */
    size_t nx;             /* size of matrix */
    size_t ny;             /* size of matrix */
    double *u_1;              /* output matrix */
    double *u_2;              /* output matrix */

    /* check for proper number of arguments */
    if(nrhs!=6) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs","Two inputs required.");
    }
    if(nlhs!=2) {
        mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nlhs","One output required.");
    }    
    
    A = mxGetPr(prhs[0]);
    B = mxGetPr(prhs[1]);
    c = mxGetPr(prhs[2]);
    d_1 = mxGetPr(prhs[3]);
    d_2 = mxGetPr(prhs[4]);
    alpha = mxGetScalar(prhs[5]);

    /* get dimensions of the input matrix */
    nx = mxGetM(prhs[1]);
    ny = mxGetN(prhs[1]);

    /* create the output matrix */
    plhs[0] = mxCreateDoubleMatrix((mwSize)nx,(mwSize)ny,mxREAL);
    plhs[1] = mxCreateDoubleMatrix((mwSize)nx,(mwSize)ny,mxREAL);

    /* get a pointer to the real data in the output matrix */
    d_1 = mxGetPr(plhs[0]);
    d_2 = mxGetPr(plhs[1]);
    
    mwSize i;
    double temp;
    double invalpha;
    double d1;
    double d2;
    double lambda1;
    double lambda2;
    double x;
    double y;
    invalpha = 1/alpha;
    
    for (i=0; i<nx*ny; i++) {
       if (A[i]==0 & B[i]==0){
           u_1[i] = d_1[i];
           u_2[i] = d_2[i];
       }
       else if (A[i]!=0 & B[i]==0){
           temp = SoftShrink(-c[i] / A[i] + d_1[i], invalpha * fabs(A[i]));
           u_1[i] = temp + c[i] / A[i];
           u_2[i] = d_2[i];        
       }
       else if (A[i]==0 & B[i]!=0){
           u_1[i] = d_1[i];
           temp = SoftShrink(-c[i] / B[i] + d_2[i], invalpha * fabs(B[i]));
           u_2[i] = temp + c[i] / B[i];
       }
       else {
           d1 = A[i] * d_1[i];
           d2 = B[i] * d_2[i] - c[i];
           lambda1 = invalpha * A[i] * A[i];
           lambda2 = invalpha * B[i] * B[i];
            
            if ((d1+d2-lambda1-lambda2) > 0) {
                x[i] = d1 - lambda1;
                y[i] = d2 - lambda2;
            } else if ((d1+d2+lambda1+lambda2) < 0) {
                x[i] = d1 + lambda1;
                y[i] = d2 + lambda2;                
            } else {
                x[i] = (d1 * lambda2 - d2 * lambda1) / (lambda1+lambda2);
                y[i] = (d2 * lambda1 - d1 * lambda2) / (lambda1+lambda2);
                u_1[i] = x / A[i];
                u_2[i] = (y + c[i]) / B[i];                
            }
       }
    }    
}
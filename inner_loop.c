//
//  code.c
//  CLoop
//
//  Created by Mathieu Zerter on 2015-02-11.
//
//


#include "mex.h"

/* The gateway function */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* variable declarations here */
    
    double *ptr;
    
    double *inMatrix1;
    double *inMatrix2;
    double *inMatrix3;
    double *inMatrix4;
	double *wvec;
    
    int ncols1, ncols2, ncols3, ncols4, nrows1, nrows2, nrows3, nrows4;
    
    ncols1 = mxGetM(prhs[0]);
    nrows1 = mxGetN(prhs[0]);
    
    ncols2 = mxGetM(prhs[1]);
    nrows2 = mxGetN(prhs[1]);
    
    ncols3 = mxGetM(prhs[2]);
    nrows3 = mxGetN(prhs[2]);
    
    ncols4 = mxGetM(prhs[3]);
    nrows4 = mxGetN(prhs[3]);
    
    plhs[0] = mxCreateDoubleMatrix(1,ncols1*ncols2*ncols3*ncols4,mxREAL);
    
    inMatrix1 = mxGetPr(prhs[0]);
    inMatrix2 = mxGetPr(prhs[1]);
    inMatrix3 = mxGetPr(prhs[2]);
    inMatrix4 = mxGetPr(prhs[3]);
    
    wvec = mxGetPr(prhs[4]);
    
    ptr = mxGetPr(plhs[0]);
    
		for (int m = 0; m < nrows1; m++)
		{

			double wv = wvec[m];
            
			int mn1 = m*ncols1;
			int mn2 = m*ncols2;
			int mn3 = m*ncols3;
			int mn4 = m*ncols4;

			for (int i = 0; i < ncols1; i++)
			{
				double arri = inMatrix1[i + mn1];

				for (int j = 0; j < ncols2; j++)
				{

					double arrj = inMatrix2[j + mn2];

					for (int k = 0; k < ncols3; k++)
					{

						double arrk = inMatrix3[k + mn3];

						for (int l = 0; l < ncols4; l++)
						{
							int index = i*ncols4*ncols2*ncols3 + j 
									* ncols4 * ncols3 + k*ncols4 + l;

							ptr[index] += arri *arrj
								* arrk * inMatrix4[l + mn4] * wv;

						}
					}



				}
			}
            
		}
    
    
        
}

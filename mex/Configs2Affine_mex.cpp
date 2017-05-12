#include "mex.h"
#include <math.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[], 
				 int nrhs, const mxArray *prhs[])
{

	int numConfigs;
	//double a11, a12, a13, a21, a22, a23;

	// input variables
	double *configs;
	double tx, ty, r2, sx, sy, r1;
	int sourceH,sourceW;
	int targetH,targetW;
	int r1x, r1y, r2x, r2y;


	// output variables
	double *affines;
	int *insiders;


	/* Find the dimensions of the data */
	numConfigs = mxGetN(prhs[0]);


	/* Create an mxArray for the output data */
	plhs[0] = mxCreateDoubleMatrix(6, numConfigs, mxREAL );
	plhs[1] = mxCreateNumericMatrix (1, numConfigs, mxINT32_CLASS, mxREAL);


	/* Retrieve the input data */

	configs = mxGetPr(prhs[0]);
	sourceH = mxGetScalar(prhs[1]);
	sourceW = mxGetScalar(prhs[2]);
	targetH = mxGetScalar(prhs[3]);
	targetW = mxGetScalar(prhs[4]);
	r1x = mxGetScalar(prhs[5]);
	r1y = mxGetScalar(prhs[6]);
	r2x = mxGetScalar(prhs[7]);
	r2y = mxGetScalar(prhs[8]);
 

	/* Create a pointer to the output data */
	affines = mxGetPr(plhs[0]);
	insiders = (int*)mxGetPr(plhs[1]);


	// MAIN LOOP
	for (int i = 0 ; i < numConfigs ; i++)
	{
		//if (i==48195)
		//	mexPrintf("MAIN LOOP\n");
// 		if (i%100000==0)
// 			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);

		tx = configs[6*i];
		ty = configs[6*i+1];
		r2 = configs[6*i+2];
		sx = configs[6*i+3];
		sy = configs[6*i+4];
		r1 = configs[6*i+5];

		//syms tx ty r2 sx sy r1;
		// A =     [ sx*cos(r1)*cos(r2) - sy*sin(r1)*sin(r2), - sx*cos(r1)*sin(r2) - sy*cos(r2)*sin(r1), tx]
		//         [ sx*cos(r2)*sin(r1) + sy*cos(r1)*sin(r2),   sy*cos(r1)*cos(r2) - sx*sin(r1)*sin(r2), ty]
		//         [                                       0,                                         0,  1]
		//        matrixConfigs(i,:) = [A(1,1),A(1,2),A(1,3),A(2,1),A(2,2),A(2,3)];

		affines[6*i]   = sx*cos(r1)*cos(r2) - sy*sin(r1)*sin(r2);
		affines[6*i+1] = - sx*cos(r1)*sin(r2) - sy*cos(r2)*sin(r1);
		affines[6*i+2] = tx;
		affines[6*i+3] = sx*cos(r2)*sin(r1) + sy*cos(r1)*sin(r2);
		affines[6*i+4] = sy*cos(r1)*cos(r2) - sx*sin(r1)*sin(r2);
		affines[6*i+5] = ty;

		//mexPrintf("%.3f, %.3f, %.3f, %.3f, %.3f, %.3f\n",a11,a12,a13,a21,a22,a23);
                
		double c1x = affines[6*i]  *(1-(r1x+1)) + affines[6*i+1]*(1-(r1y+1)) + (r2x+1) + tx;
		double c1y = affines[6*i+3]*(1-(r1x+1)) + affines[6*i+4]*(1-(r1y+1)) + (r2y+1) + ty;
		double c2x = affines[6*i]  *(sourceW-(r1x+1)) + affines[6*i+1]*(1-(r1y+1)) + (r2x+1) + tx;
		double c2y = affines[6*i+3]*(sourceW-(r1x+1)) + affines[6*i+4]*(1-(r1y+1)) + (r2y+1) + ty;
		double c3x = affines[6*i]  *(sourceW-(r1x+1)) + affines[6*i+1]*(sourceH-(r1y+1)) + (r2x+1) + tx;
		double c3y = affines[6*i+3]*(sourceW-(r1x+1)) + affines[6*i+4]*(sourceH-(r1y+1)) + (r2y+1) + ty;
		double c4x = affines[6*i]  *(1-(r1x+1)) + affines[6*i+1]*(sourceH-(r1y+1)) + (r2x+1) + tx;
		double c4y = affines[6*i+3]*(1-(r1x+1)) + affines[6*i+4]*(sourceH-(r1y+1)) + (r2y+1) + ty;

                // allow to exceed boundary by at most 10
		int k = int((c1x>-10)&&(c1x<targetW+10)&&(c1y>-10)&&(c1y<targetH+10)&&
                    (c2x>-10)&&(c2x<targetW+10)&&(c2y>-10)&&(c2y<targetH+10)&&
                    (c3x>-10)&&(c3x<targetW+10)&&(c3y>-10)&&(c3y<targetH+10)&&
                    (c4x>-10)&&(c4x<targetW+10)&&(c4y>-10)&&(c4y<targetH+10));
		//int kk = (5==5);
		//int kkk = (5==4);
		//if (i%100==0)
		//	mexPrintf("%d,%d,%d\n",k,kk,kkk);
		insiders[i] = k;
	}

}
//
//%
//cornersX = [1 w1 w1 1];
//cornersY = [1 1 h1 h1];
//cornersA = round(a2x2*[cornersX-(r1x+1);cornersY-(r1y+1)]);
//cornerAxs = round(cornersA(1,:) + (r2x+1)  + a(1,3));
//cornerAys = round(cornersA(2,:) + (r2y+1)  + a(2,3));

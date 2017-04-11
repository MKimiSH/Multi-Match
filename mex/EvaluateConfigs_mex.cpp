#include "mex.h"
#include <math.h>
#include <stdio.h>
#include <string.h>

#ifndef max
#define max( a, b ) ( ((a) > (b)) ? (a) : (b) )
#endif

#ifndef min
#define min( a, b ) ( ((a) < (b)) ? (a) : (b) )
#endif

void mexFunction(int nlhs, mxArray *plhs[],
        int nrhs, const mxArray *prhs[]) {
    // parameters
    int h1, w1, h2, w2;
    int r1x, r1y, r2x, r2y;
    
    int numPoints, numConfigs;
    double a11, a12, a13, a21, a22, a23;
    
    // input variables
    double *img1, *img2;
    double *configs;
    int *xs, *ys;
    
    // Should we use photometrics?
    int usePhotometrics;
    
    // helper variables
    int *xs_centered, *ys_centered;
    double *valsI1;
    int targetPoint_x, targetPoint_y;
    int targetInd;
    double score;
    int maxInd2img2;
    
    double sigX;
    double sigY;
    double meanX;
    double meanY;

    
    // output variables
    double *distances;
    
    
    //// checking inputs
    //if (nrhs != 5)
    //	mexErrMsgTxt("The number of input arguments must be 5.");
    //if (nlhs != 1)
    //	mexErrMsgTxt("The number of output arguments must be 5.");
    
    
    // inputs:     dimensions
    // ------	    -----------
    // image I1     h1 x w1
    // image I2     h2 x w2
    // configs      numConfigs x 6
    // xs           1 x numPoints
    // ys1			1 x numPoints
    
    // outputs:     dimensions
    // -------	    -----------
    // distances    1 x numConfigs
    
    
    /* Find the dimensions of the data */
    h1 = mxGetN(prhs[0]);
    w1 = mxGetM(prhs[0]);
    h2 = mxGetN(prhs[1]);
    w2 = mxGetM(prhs[1]);
    numConfigs = mxGetN(prhs[2]);
    numPoints = mxGetN(prhs[3]);
    
    r1x = 0.5*(w1-1);
    r1y = 0.5*(h1-1);
    r2x = 0.5*(w2-1);
    r2y = 0.5*(h2-1);
    
    maxInd2img2 = h2*w2 - 1;
    
    
    /* Create an mxArray for the output data */
    plhs[0] = mxCreateDoubleMatrix(1, numConfigs, mxREAL );
    
    /* Create an mxArrays for temporary data */
    xs_centered = (int *)malloc(numPoints*sizeof(int));
    ys_centered = (int *)malloc(numPoints*sizeof(int));
    valsI1 = (double *)malloc(numPoints*sizeof(double));
    
    /* Store target pixel locations - x and y */
    double* xs_target = (double *)malloc(numPoints*sizeof(double));
    double* ys_target = (double *)malloc(numPoints*sizeof(double));
    
    
    /* Retrieve the input data */
    img1 = mxGetPr(prhs[0]);
    double* tmp_img2 = mxGetPr(prhs[1]);
    configs = mxGetPr(prhs[2]);
    xs = (int*)mxGetPr(prhs[3]);
    ys = (int*)mxGetPr(prhs[4]);
    usePhotometrics = mxGetScalar(prhs[5]);
    // img2 is of height 3*img2 (this padding is for not needing to check bounds)
    img2 = (double*)malloc(5*h2*w2*sizeof(double));
    memset(img2, 2, 5*h2*w2*sizeof(double));
    memcpy(img2+2*h2*w2, tmp_img2, h2*w2*sizeof(double));
    //
    
    int cN = mxGetN(prhs[2]);
    int cM = mxGetM(prhs[2]);
    
    
    /*Centered pointes*/
    for (int i = 0 ; i < numPoints ; i++) {
        xs_centered[i] = xs[i]-(r1x+1);
        ys_centered[i] = ys[i]-(r1y+1);
    }
    
    /*Precalculating source point indices into I1 (and the values themselves)*/
    for (int j = 0; j < numPoints ; j++) {
        //sourceInds[j] = (ys[j] - 1)*w1 + xs[j];
        valsI1[j] = img1[(ys[j] - 1)*w1 + xs[j]-1]; // -1 is for c
    }
    
    
    /* Create a pointer to the output data */
    distances = mxGetPr(plhs[0]);

// plhs[1] = mxCreateDoubleMatrix(1, numConfigs, mxREAL);
// plhs[2] = mxCreateDoubleMatrix(1, numConfigs, mxREAL);
// plhs[3] = mxCreateDoubleMatrix(1, numConfigs, mxREAL);
// plhs[4] = mxCreateDoubleMatrix(1, numConfigs, mxREAL);
// double* outMeanX = mxGetPr(plhs[1]);
// double* outMeanY = mxGetPr(plhs[2]);
// double* outSigX = mxGetPr(plhs[3]);
// double* outSigY = mxGetPr(plhs[4]);

        
    // MAIN LOOP
    
    for (int i = 0 ; i < numConfigs ; i++) {
        
//         if (i%100000==0)
//             mexPrintf("MAIN LOOP: config %d out of %d\n", i+1, numConfigs);
//		if (i==709416)
//			mexPrintf("MAIN LOOP: config %d out of %d\n",i+1,numConfigs);
        
        a11 = configs[6*i];
        a12 = configs[6*i+1];
        a13 = configs[6*i+2];
        a21 = configs[6*i+3];
        a22 = configs[6*i+4];
        a23 = configs[6*i+5];
        
        
        double score = 0;
        
        double tmp1 = (r2x+1) + a13 + 0.5;
        double tmp2 = (r2y+1) + a23 + 0.5 + 2*h2;
        double* ptrVals = valsI1;
        int* ptrXsc = xs_centered;
        int* ptrYsc = ys_centered;
        if (!usePhotometrics) {
            for (int j = 0; j < numPoints ; j++) {
                targetPoint_x = int(a11*(*ptrXsc) + a12*(*ptrYsc) + tmp1); // includes rounding
                targetPoint_y = int(a21*(*ptrXsc) + a22*(*ptrYsc) + tmp2); // includes rounding
                targetInd = (targetPoint_y - 1)*w2 + targetPoint_x - 1; // -1 is for c
                score += (fabs((*ptrVals) - img2[targetInd]));
                
                ptrVals++;
                ptrXsc++;
                ptrYsc++;
            }
        } else {
            double sumXiYi =0; double sumXi =0; double sumYi =0;
            double sumXiSqrd =0; double sumYiSqrd =0; double Xi, Yi;
            
            for (int j = 0; j < numPoints ; j++) {
                targetPoint_x = int(a11*(*ptrXsc) + a12*(*ptrYsc) + tmp1); // includes rounding
                targetPoint_y = int(a21*(*ptrXsc) + a22*(*ptrYsc) + tmp2); // includes rounding
                targetInd = (targetPoint_y - 1)*w2 + targetPoint_x - 1; // -1 is for c
                Xi = (*ptrVals) ;
                // Yi = img2[targetInd];
                xs_target[j] = Xi;
                ys_target[j] = Yi;
                
//                                 sumXiYi += Xi*Yi;
                sumXi += Xi;
                sumYi += Yi;
                sumXiSqrd += (Xi*Xi);
                sumYiSqrd += (Yi*Yi);
                ptrVals++;
                ptrXsc++;
                ptrYsc++;
            }
            // old score based on correlation
//                         score = numPoints*sumXiYi-sumXi*sumYi;
//                         score /= (0.00001 + sqrt(numPoints*sumXiSqrd - sumXi*sumXi) * sqrt(numPoints*sumYiSqrd - sumYi*sumYi) );
//                         score  = (1-score)/2;
// 			score *= numPoints;
            // new score, based on normalizing mean and std of each of the signals
            double epsilon = 0.0000001;
            sigX = sqrt((sumXiSqrd-(sumXi*sumXi)/numPoints)/numPoints) + epsilon;
            sigY = sqrt((sumYiSqrd-(sumYi*sumYi)/numPoints)/numPoints) + epsilon;
            meanX = sumXi/numPoints;
            meanY = sumYi/numPoints;
            double sigXoversigY = sigX/sigY;
            
            // Variable that stores a sum used repeatadly in the computation: -meanX+sigXoversigY*meanY
            double faster = -meanX+sigXoversigY*meanY;
            
            for (int j = 0; j < numPoints ; j++) {
                //                       score += fabs( (xs_target[j]-meanX)/sigX - (ys_target[j]-meanY)/sigY);
//                                score += fabs( (xs_target[j]-meanX) - sigX*(ys_target[j]-meanY)/sigY);
//                                score += fabs( (xs_target[j]-meanX) - sigXoversigY*(ys_target[j]-meanY));
                score += fabs( xs_target[j]- sigXoversigY*ys_target[j] + faster);
            }
            
        }
        /*intermediate*/
        //for (int j = 0; j < numPoints ; j++)
        //{
        //	//targetPoint_x = int(a11*xs_centered[j] + a12*ys_centered[j] + (r2x+1) + a13 + 0.5); // includes rounding
        //	//targetPoint_y = int(int(0.5*h2) + a21*xs_centered[j] + a22*ys_centered[j] + (r2y+1) + a23 + 0.5); // includes rounding
        //	targetPoint_x = int(a11*xs_centered[j] + a12*ys_centered[j] + tmp1); // includes rounding
        //	targetPoint_y = int(int(0.5*h2) + a21*xs_centered[j] + a22*ys_centered[j] + tmp2); // includes rounding
        //	targetInd = (targetPoint_y - 1)*w2 + targetPoint_x - 1; // -1 is for c
        //	//targetInd = max(0,min(targetInd,maxInd2img2));
        //	//targetInd = max(0,targetInd);
        //	score += (fabs(valsI1[j] - img2[targetInd]));
        //}
        /*original*/
        //for (int j = 0; j < numPoints ; j++)
        //{
        //	targetPoint_x = int(a11*xs_centered[j] + a12*ys_centered[j] + (r2x+1) + a13 + 0.5); // includes rounding
        //	targetPoint_y = int(a21*xs_centered[j] + a22*ys_centered[j] + (r2y+1) + a23 + 0.5); // includes rounding
        //	targetInd = (targetPoint_y - 1)*w2 + targetPoint_x - 1; // -1 is for c
        //	targetInd = max(0,min(targetInd,maxInd2img2));
        //	score += (fabs(valsI1[j] - img2[targetInd]));
        //}
        
        distances[i] = score/numPoints;

// outMeanX[i] = meanX;
// outMeanY[i] = meanY;
// outSigX[i] = sigX;
// outSigY[i] = sigY;
    }
    
    /* Free the allocated arrays */
    free(xs_centered);
    free(ys_centered);
    free(valsI1);
    free(img2);
    
}
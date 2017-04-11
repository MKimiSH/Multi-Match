//function configs = CreateList(ntx_steps,nty_steps,nr_steps,NR2_steps,ns_steps,tx_steps,ty_steps,r_steps,s_steps,gridSize)
//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
// configs = zeros(gridSize,6);
// gridInd = 0;
// for tx_ind = 1 : ntx_steps 
//     tx = tx_steps(tx_ind);
// %     fprintf('step %d out of %d\n',tx_ind,nt_steps);
//     for ty_ind = 1 : nty_steps 
//         ty = ty_steps(ty_ind);
//         for r1_ind = 1 : nr_steps
//             r1 = r_steps(r1_ind);
//             for r2_ind =  1 : NR2_steps
//                 r2 = r_steps(r2_ind);
//                 for sx_ind =  1 : ns_steps
//                     sx = s_steps(sx_ind);
//                     for sy_ind =  1 : ns_steps 
//                         sy = s_steps(sy_ind);
//                         
//                         gridInd = gridInd + 1;
//                         
//                         configs(gridInd,:) = [tx,ty,r2,sx,sy,r1];
//                     end
//                 end
//             end
//         end
//     end
// end


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

	/* Retrieve the input data */
	int ntx_steps = mxGetScalar(prhs[0]);
	int nty_steps = mxGetScalar(prhs[1]);
	int nr_steps = mxGetScalar(prhs[2]);
	int NR2_steps = mxGetScalar(prhs[3]);
	int ns_steps = mxGetScalar(prhs[4]);
	double *tx_steps = mxGetPr(prhs[5]);
	double *ty_steps = mxGetPr(prhs[6]);
	double *r_steps = mxGetPr(prhs[7]);
	double *s_steps = mxGetPr(prhs[8]);
	int numConfigs = mxGetScalar(prhs[9]);


	/* Create an mxArray for the output data */
	plhs[0] = mxCreateDoubleMatrix(6, numConfigs, mxREAL );

	/* Create a pointer to the output data */
	double *configs = mxGetPr(plhs[0]);

	double tx, ty, r2, sx, sy, r1;

	// MAIN LOOP
	int gridInd = 0;
	for (int tx_ind = 0 ; tx_ind < ntx_steps ; tx_ind++)
	{
		tx = tx_steps[tx_ind];
		for (int ty_ind = 0 ; ty_ind < nty_steps ; ty_ind++)
		{
			ty = ty_steps[ty_ind];
			for (int r1_ind = 0 ; r1_ind < nr_steps ; r1_ind++)
			{
				r1 = r_steps[r1_ind];
				for (int r2_ind = 0 ; r2_ind < NR2_steps ; r2_ind++)
				{
					r2 = r_steps[r2_ind];
					for (int sx_ind = 0 ; sx_ind < ns_steps ; sx_ind++)
					{
						sx = s_steps[sx_ind];
						for (int sy_ind = 0 ; sy_ind < ns_steps ; sy_ind++)
						{
							sy = s_steps[sy_ind];
							// inside stuff                        
							configs[gridInd++] = tx;
							configs[gridInd++] = ty;
							configs[gridInd++] = r2;
							configs[gridInd++] = sx;
							configs[gridInd++] = sy;
							configs[gridInd++] = r1;
						}
					}
				}
			}
		}
	}

}
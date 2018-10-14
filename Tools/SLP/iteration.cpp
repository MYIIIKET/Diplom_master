#include "mex.h"
#include "matrix.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

/*
The function iteration propagate the label through the graph.
Corresponding to line 3 to line 17 in Algorithm 1
	%     De-Ming Liang and Yu-Feng Li. Lightweight Label Propagation for
	%     Large-Scale Network Data. In: Proceedings of the 27th International
	%     Joint Conference on Artificial Intelligence (IJCAI'18).
@W:			sparse matrix
@f:			predicting vector, the labeled ones are placed at the beginning end
@lnonzeros:	nonzeros' indexes in cells
@stepSize:	step size
@steps:		running epoches
*/



void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
	if (nrhs != 6) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
			"Six inputs required.");
	}

	srand((unsigned)time(NULL));
	double *f, *ftmp, *W;
	mwIndex *Ir, *Jc;
	const mxArray *lnonzeros;
	double stepSize, stepSize1, steps, l;
	size_t Wcol, Wrow, frow, fcol;

	ftmp = mxGetPr(prhs[1]);
	frow = mxGetM(prhs[1]);
	fcol = mxGetN(prhs[1]);
	f = (double *)mxCalloc(frow, sizeof(double));
	for (int i = 0; i < frow; i++) {
		f[i] = ftmp[i];
	}
	lnonzeros = prhs[2];
	stepSize1 = stepSize = mxGetScalar(prhs[3]);
	steps = mxGetScalar(prhs[4]);
	l = mxGetScalar(prhs[5]);
	double var = 1.0;
	W = mxGetPr(prhs[0]);
	
	if (mxIsSparse(prhs[0])) { // sparse matrix
		Ir = mxGetIr(prhs[0]);
		Jc = mxGetJc(prhs[0]);
		Wcol = mxGetNzmax(prhs[0]);
		for (mwIndex it = 0; it < l; it++) { // pre-propagate
			//mwIndex it = rand() % frow + 1;

			mxArray *cell_ele = mxGetCell(lnonzeros, it);
			if (cell_ele != NULL) {
				//double *indexs = mxGetPr(cell_ele);  // something wrong
				mwIndex *indexs = (mwIndex *)mxGetData(cell_ele);

				size_t M = mxGetM(cell_ele);
				size_t indNum = mxGetN(cell_ele);

				for (mwSize j = 0; j < indNum; j++) {
					mwIndex ind = indexs[j] - 1; // matlab to C version index
					mxAssert(Jc[it] + j < Wcol, "ind exceed! %d,%d\n", ind, it);
					f[ind] = f[ind] - stepSize*W[Jc[it] + j] * (f[ind] - f[it]); 
				}
			}
		}

		mwIndex it = rand() % frow;
		for (mwSize i = 0; i < steps; i++) {
			for (mwIndex k = 0; k < frow; k++) { 
				//printf("%d,%d:\n", i, it);
				it = (it + 1) % frow;
				mxArray *cell_ele = mxGetCell(lnonzeros, it);
				if (cell_ele != NULL) {
					mwIndex *indexs = (mwIndex *)mxGetData(cell_ele);

					size_t M = mxGetM(cell_ele);
					size_t indNum = mxGetN(cell_ele);

					for (mwSize j = 0; j < indNum; j++) {
						mwIndex ind = indexs[j] - 1; // tranform matlab indexes to C indexes
						//mexPrintf("%d,%d\n", ind, it);
						mxAssert(Jc[it] + j < Wcol, "ind exceed! %d,%d\n", ind, it);
						f[ind] = f[ind] - stepSize*W[Jc[it] + j] * (f[ind] - f[it]); 
					}
				}
				stepSize = stepSize1 / sqrt(var++); // 72 line
			}
		}
	}
	else { // full matrix
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
			"The first parameter should be sparse matrix.");
	}
	mxArray *re= mxCreateDoubleMatrix(frow, fcol, mxREAL);
	mxSetPr(re, f);
	plhs[0] = re;
	nlhs = 1;

	return;
}

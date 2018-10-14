#include "mex.h"
#include "matrix.h"

/*
The function find2cells find neighborhood indexes for each node and return a cells.
Corresponding to line 2 in Algorithm 1 
	%     De-Ming Liang and Yu-Feng Li. Lightweight Label Propagation for
	%     Large-Scale Network Data. In: Proceedings of the 27th International
	%     Joint Conference on Artificial Intelligence (IJCAI'18).
@row:	row indexs
@col:	col index
@edges: num of edges
@n:		num of nodes
@l:		labels nums
*/

mxArray* main_process(double *row, double *col, mwSize edges, mwSize n, mwSize l) {
	mxArray *re = mxCreateCellMatrix(n, 1);
	mwSize count = 0, cur = 1;

	for(mwSize i = 0; i < edges; i++){
		if (col[i] != cur) {
			//mexPrintf("%d\n", cur);
			mxArray *tmp = mxCreateNumericMatrix(0, 0, mxINT64_CLASS, mxREAL);
			mwSize *dy = (mwSize *)mxMalloc(count * sizeof(mwSize));
			for (mwSize j = 0; j < count; j++) {
				dy[j] = (mwSize)row[i-count+j] + l;
			}
			//mexPrintf("\n");
			mxSetPr(tmp, (double *)dy);
			mxSetM(tmp, 1);
			mxSetN(tmp, count);
			mxSetCell(re, cur-1, tmp);

			count = 1;
			//cur++;
			cur = col[i];
		}
		else {
			count++;
		}
	}
	//mexPrintf("loop out\n");
	mxArray *tmp = mxCreateNumericMatrix(0, 0, mxINT64_CLASS, mxREAL);
	mwSize *dy = (mwSize *)mxMalloc(count * sizeof(mwSize));
	for (mwSize j = 0; j < count; j++) {
		dy[j] = (mwSize)row[edges - count + j] + l;
	}
	mxSetPr(tmp, (double *)dy);
	mxSetM(tmp, 1);
	mxSetN(tmp, count);
	mxSetCell(re, cur-1, tmp);
	return re;
}


void mexFunction(mwSize nlhs, mxArray *plhs[], mwSize nrhs, const mxArray *prhs[])
{
	if (nrhs != 5) {
		mexErrMsgIdAndTxt("MyToolbox:arrayProduct:nrhs",
			"Five inputs required.");
	}
	
	
	double *row, *col;
	mwSize edges, n, l;


	row = mxGetPr(prhs[0]);
	col = mxGetPr(prhs[1]);
	edges = mxGetScalar(prhs[2]);
	n = mxGetScalar(prhs[3]);
	l = mxGetScalar(prhs[4]);

	plhs[0] = main_process(row, col, edges, n, l);
	nlhs = 1;

	return;
}

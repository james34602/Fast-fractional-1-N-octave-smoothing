#include <stdio.h>
#include <math.h>
#include <float.h>
#include <mex.h>
#if !defined(MAX)
#define MAX(A, B) ((A) > (B) ? (A) : (B))
#endif
#if !defined(MIN)
#define MIN(A, B) ((A) < (B) ? (A) : (B))
#endif
#include "interpolation2.h"
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
	// Check for proper number of arguments
	if (nrhs != 2)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:invalidNumInputs", "Please specify struct and spectrum data.");
	else if (nlhs > 1)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:maxlhs", "Too many output arguments.");
	char *octaveSmooth = (char*)mxGetPr(prhs[0]);
	size_t dm1 = mxGetM(prhs[0]);
	size_t dn1 = mxGetN(prhs[0]);
	if (MIN(dm1, dn1) != 1)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:maxlhs", "Input struct must be a vector.");
	double *spectrum = (double*)mxGetPr(prhs[1]);
	size_t dm2 = mxGetM(prhs[1]);
	size_t dn2 = mxGetN(prhs[1]);
	size_t m2 = MAX(dm2, dn2);
	if (MIN(dm2, dn2) != 1)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:maxlhs", "Input spectrum must be a vector.");
	// Create a matrix for the return argument
	if (m2 == dm2)
		plhs[0] = mxCreateNumericMatrix(m2, 1, mxDOUBLE_CLASS, mxREAL);
	else
		plhs[0] = mxCreateNumericMatrix(1, m2, mxDOUBLE_CLASS, mxREAL);
	// Assign pointers to the various parameters
	double *Y = (double*)mxGetPr(plhs[0]);
	// Do the actual computations in a subroutine
	unsigned int aryLen = *((unsigned int*)octaveSmooth);
    if (aryLen < 1)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:maxlhs", "Something very wrong with struct vector");
	size_t dataStructSize = sizeof(ierper) + sizeof(double) * aryLen * 2 + sizeof(double) * 3 + (sizeof(double) * (aryLen - 4) * 4) + sizeof(double) * 3;
	double *weighting = (double*)(octaveSmooth + sizeof(unsigned int));
	unsigned int *efIdx = (unsigned int*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double));
	ierper *interpolator = (ierper*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double) + (aryLen + 1) * sizeof(unsigned int));
	unsigned int recNout = *((unsigned int*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double) + (aryLen + 1) * sizeof(unsigned int) + dataStructSize));
	unsigned int skBin = *((unsigned int*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double) + (aryLen + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int) + sizeof(double) * recNout + sizeof(double) * aryLen));
    if (recNout < 1)
		mexErrMsgIdAndTxt("MATLAB:james_smoothOctave:maxlhs", "Something very wrong with struct vector");
	double *freqArray = (double*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double) + (aryLen + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int));
	double *Px_oct = (double*)(octaveSmooth + sizeof(unsigned int) + aryLen * sizeof(double) + (aryLen + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int) + sizeof(double) * recNout);
    unsigned int i;
	for (i = 0; i < aryLen; i++)
	{
		double mean = 0.0;
		for (unsigned int j = efIdx[i]; j <= efIdx[i + 1]; j++)
			mean += spectrum[j];
		Px_oct[i] = mean * weighting[i];
	}
	makimaUpdate(interpolator, Px_oct, 1, 1);
    if (m2 != recNout)
        printf("Input spectrum size doesn't match the one you provide in the initialization stage\n");
	for (i = 0; i < skBin; i++)
        Y[i] = spectrum[i];
	for (i = skBin; i < recNout; i++)
		Y[i] = interpolator->interpolate(interpolator, freqArray[i], Px_oct);
	return;
}
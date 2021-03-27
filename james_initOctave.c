#include <stdio.h>
#include <stdlib.h>
#include <string.h>
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
unsigned int auditoryInterpIndex(const double nocts, double *freq, const unsigned int len, const unsigned int fcLen, double *fc, double *edgeFreq, unsigned int *edgeFreqIdx)
{
	unsigned int i, j;
	// octave smoothing
	double nocts2x = nocts * 2.0;
	double f1 = 1.0;
	for (i = 0; i < fcLen; i++)
	{
		f1 = f1 * pow(10.0, 3.0 / (10.0 * nocts2x));
		fc[i] = f1;
	}
	for (i = 0; i < fcLen; i++)
		edgeFreq[i] = fc[i] * pow(10, 3.0 / (20.0 * nocts2x));
	for (i = 0; i < fcLen; i++)
	{
		int fe_p = -1, fe_m = -1, fe_0 = -1;
		for (j = 0; j < len; j++)
		{
			if (freq[j] > edgeFreq[i])
			{
				fe_p = j;
				break;
			}
		}
		for (j = 0; j < len; j++)
		{
			if (freq[j] >= edgeFreq[i])
			{
				break;
			}
			else
				fe_m = j;
		}
		for (j = 0; j < len; j++)
		{
			if (freq[j] == edgeFreq[i])
			{
				fe_0 = j;
				break;
			}
		}
		if (fe_0 > -1)
			edgeFreq[i] = fe_0 + 1.0;
		else
		{
			double p = (fe_p + 1.0) - edgeFreq[i];
			double m = edgeFreq[i] - (fe_m + 1.0);
			if (p < m && fe_p > -1)
				edgeFreq[i] = fe_p + 1.0;
			else
				edgeFreq[i] = fe_m + 1.0;
		}
	}
	for (i = 0; i < fcLen; i++)
		edgeFreqIdx[i] = ((unsigned int)edgeFreq[i]) - 1; // Convert to C indexing
	for (i = 1; i < fcLen; i++)
		fc[i - 1] = fc[i]; // Shift
	fc[fcLen - 1] = -1.0; // Invalid index notation
	unsigned int skipBin = 0;
	for (i = 0; i < fcLen; i++)
	{
		if (fc[0] < freq[i])
		{
			skipBin = i;
			break;
		}
	}
	return skipBin;
}
unsigned int getAuditoryBandLen(const double nocts, double *freq, const unsigned int len)
{
	if (nocts < 0.15)
		return 0;
	// octave smoothing
	double nocts2x = nocts * 2.0;
	// octave center frequencies
	double f1 = 1.0;
	unsigned int i = 0;
	while (f1 < freq[len - 1])
	{
		f1 = f1 * pow(10.0, 3.0 / (10.0 * nocts2x));
		i++;
	}
	return i;
}
void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    // Check for proper number of arguments
    if (nrhs != 2)
        mexErrMsgIdAndTxt("MATLAB:james_initOctave:invalidNumInputs", "Please specify 1/N-Octave and frequency vector");
    else if (nlhs > 1)
        mexErrMsgIdAndTxt("MATLAB:james_initOctave:maxlhs", "Too many output arguments.");
	size_t dm1 = mxGetM(prhs[0]);
	size_t dn1 = mxGetN(prhs[0]);
	if (dm1 != dn1 || dm1 > 1 || dn1 > 1)
		mexErrMsgIdAndTxt("MATLAB:james_initOctave:maxlhs", "N must be a double");
	double *X = (double*)mxGetPr(prhs[0]);
    double Noct = *X;
    if (Noct < 0.15)
		mexErrMsgIdAndTxt("MATLAB:james_initOctave:maxlhs", "Invalid 1/N octave");
	size_t dm2 = mxGetM(prhs[1]);
	size_t dn2 = mxGetN(prhs[1]);
	size_t Nout = MAX(dm2, dn2);
	if (MIN(dm2, dn2) != 1)
		mexErrMsgIdAndTxt("MATLAB:james_initOctave:maxlhs", "Frequency array must be a vector.");
	double *freq = (double*)mxGetPr(prhs[1]);
    // Do the actual computations in a subroutine
	unsigned int fcLen = getAuditoryBandLen(Noct, freq, Nout);
	double *fc_oct = (double*)malloc(fcLen * sizeof(double));
	double *edgeFreq = (double*)malloc(fcLen * sizeof(double));
	unsigned int *edgeFreqIdx = (unsigned int*)malloc(fcLen * sizeof(unsigned int));
	unsigned int skipBin = 0;
	if (fcLen)
		skipBin = auditoryInterpIndex(Noct, freq, Nout, fcLen, fc_oct, edgeFreq, edgeFreqIdx);
    else
		mexErrMsgIdAndTxt("MATLAB:james_initOctave:maxlhs", "Possibly invalid sample rate.");
	free(edgeFreq);
	unsigned int m1 = fcLen - 1;
	double *multiplicationPrecompute = (double*)malloc(m1 * sizeof(double));
    unsigned int i;
	for (i = 0; i < m1; i++)
		multiplicationPrecompute[i] = 1.0 / (edgeFreqIdx[i + 1] - edgeFreqIdx[i] + 1.0);
	size_t dataStructSize = sizeof(ierper) + sizeof(double) * m1 * 2 + sizeof(double) * 3 + (sizeof(double) * (m1 - 4) * 4) + sizeof(double) * 3;
	size_t virtualStructSize = sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int) + sizeof(double) * Nout + sizeof(double) * m1 + sizeof(unsigned int);
    // Create a matrix for the return argument
    plhs[0] = mxCreateNumericMatrix(1, virtualStructSize, mxUINT8_CLASS, mxREAL);
    // Assign pointers to the various parameters
	char *octaveSmooth = (char*)mxGetPr(plhs[0]);
	ierper *yp = (ierper*)(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int));
	makimaPC(yp, fc_oct, m1);
	free(fc_oct);
	unsigned int *val = (unsigned int*)octaveSmooth;
	*val = m1;
	memcpy(octaveSmooth + sizeof(unsigned int), multiplicationPrecompute, m1 * sizeof(double));
	memcpy(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double), edgeFreqIdx, (m1 + 1) * sizeof(unsigned int));
	memcpy(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int), yp, dataStructSize);
	val = (unsigned int*)(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int) + dataStructSize);
	*val = Nout;
	memcpy(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int), freq, Nout * sizeof(double));
	val = (unsigned int*)(octaveSmooth + sizeof(unsigned int) + m1 * sizeof(double) + (m1 + 1) * sizeof(unsigned int) + dataStructSize + sizeof(unsigned int) + sizeof(double) * Nout + sizeof(double) * m1);
	*val = skipBin;
	free(edgeFreqIdx);
	free(multiplicationPrecompute);
    return;
}
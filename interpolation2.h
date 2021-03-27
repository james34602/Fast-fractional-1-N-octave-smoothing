typedef struct mmLerp
{
	double(*interpolate)(struct mmLerp*, double, double *);
	int n;
} ierper;
extern void makimaPC(ierper *intp, double * x, int n);
extern void makimaUpdate(ierper *intp, double * y, int left_endpoint_derivative, int right_endpoint_derivative);
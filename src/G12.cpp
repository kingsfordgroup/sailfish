#include "nnls.h"

///Compute orthogonal rotation matrix 
///Compute matrix (c, s) so that (c, s)(A) = (SQRT(A ** 2 + B ** 2))
///				  (-s, c) 		(-s, c)(B) = (	0		     )   
///Compute SIG = SQRT(A ** 2 + B ** 2)
///SIG is computed last to allow for the possibility that  SIG may be in the same location as A or B.
void G1(double* A, double* B, double* COS, double* SIN, double* SIG)
{
	const int Zero = 0;
	const int One = 1;
	double XR, YR;
	
	if (dAbs(*A) > dAbs(*B))
	{
		XR = *B / *A;
		YR = (double)sqrt(One + pow(XR, 2));
		*COS = dSign(One / YR, *A);
		*SIN = (*COS) * XR;
		*SIG = dAbs(*A) * YR;
		return;
	}
	else
	{
		if (*B != 0)
		{
			XR = *A / *B;
			YR = (double)sqrt(One + pow(XR, 2));	
			*SIN = dSign(One / YR, *B);
			*COS = (*SIN) * XR;
			*SIG = dAbs(*B) * YR;
			return;
		}
		else
		{
			*SIG = (double)Zero;
			*COS = (double)Zero;
			*SIN = (double)One;
			return;
		}
	}
}

///Apply the rotation computed by G1 to (X, Y)
void G2(double* COS, double* SIN, double* X, double* Y)
{
	double XR;
	
	XR = (*COS) * (*X) + (*SIN) * (*Y);
	*Y = -(*SIN) * (*X) + (*COS) * (*Y);
	*X = XR;
	return;
}

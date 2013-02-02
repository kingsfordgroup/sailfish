#include "nnls.h"

///Construction and/or application of a single Householder transformation.  Q = I + U *(U * * T) / B
///Mode				= 1 or 2 To select Algorithm H1 or H2
///LPivot 			is the index of the pivot element
///L1, M 			If L1 <= M, the transformation will be constructed to zero elements indexed from L1 through M
///					If L1 > M, the subroutine does an identity transformation
///U() , IUE, Up	On entry to H1 U() containes the pivot vector. IUE is the storage increment between elements.
///					On exit from H1 U() and Up contain quantities defining the vector U of the householder transformation.
///					On entry to H2 U() and Up should contain quantities previously computed by H1. These will not be modified by H2.
///C() 				On entry to H1 or H2 C() contains a matrix which will be regarded as a set of vectors to which the householder transformation is to be applied.
///					On exit C() contains the set of transformated vectors
///N				The column number of the matrix
///ICE				Storage increacement between elements of vectors in C().
///ICV				Storage increment between vectors in C().
///NCV				Number of vectors in C() to be tranformed. If NCV <= 0, no operation will be done on C().

void H12(int Mode, int LPivot, int L1, int M, double *U, int IUE, double* Up, double *C, int ICE, int ICV, int NCV)
{
	const int One = 1;
	double CL, CLInverse;
	double SM, B;
	int i, j;
	int	Incr, Index2, Index3, Index4;
	
	if ((LPivot <= 0) || (LPivot >= L1) || (L1 > M))
		return;
		
	CL = dAbs(U[LPivot]);
	if (Mode == 1)
	{
		for (j = L1; j <= M; j++)
			CL = dMax(dAbs(U[j]), CL);
			
		if (CL <= 0)
			return;
		
		CLInverse = One / CL;
		SM = (double)pow(U[LPivot] * CLInverse, 2);
		for (j = L1; j <= M; j++)
			SM = SM + (double)pow(U[j] * CLInverse , 2);
			
		CL = CL * ((double)sqrt(SM));
		if (U[LPivot] > 0)
			CL = -CL;
		(*Up) = U[LPivot] - CL;
		U[LPivot] = CL;	
	}
	else if ((Mode == 2) && (CL <= 0))
		return;
		
	if (NCV <= 0)
		return;
		
	B = (*Up) * U[LPivot];
	
	if (B >= 0)
		return;
	B = One / B;	
	
	Index2 = 1 - ICV + ICE * (LPivot - 1);
	Incr =  ICE * (L1 - LPivot);
	
	for (j = 1; j <= NCV; j++)
	{		
		Index2 = Index2 + ICV;
		Index3 = Index2 + Incr;
		Index4 = Index3;
		
		SM = C[Index2] * (*Up);
		for (i = L1; i <= M; i++)
		{
			SM = SM + C[Index3] * U[i];
			Index3 = Index3 + ICE;
		}
				
		if (SM != 0)
		{
			SM = SM * B;
			C[Index2] = C[Index2] + SM * (*Up);
			for (i = L1; i <= M; i++)
			{
				C[Index4] = C[Index4] + SM * U[i];					
				Index4 = Index4 + ICE;
			}				
		}
	}
}

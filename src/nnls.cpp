#include "nnls.h"

///              ******Nonnegative Least squares******
///Given An M by N matrix A, and an M-vector B. Compute an N-vector X, which solves the least squares problem.
///A*X=B subject to X >= 0
///
///A(), MDA, M, N	MDA is the first dimensioning parameter for the array A().
///				On entry A() contains the M by N matrix A. On exit A() contains the product matrix Q*A. 
///				Where Q is an M by M orthogonal matrix generated implicitly by this subroutine.
///B()			On entry B() contains the M-vecotr B. On exit B() contains  Q*B.
///X()			On entry X() need not be initialized. On exit X() will contain the solution vector.
///RNORM		On exit RNORM contains the euclidean norm of the residual vector.
///W()			An N-array of working space. On exit W() will contain the dual solution vector. 	
///				W() will satisfy W(I) = 0 for all I in set P and W(I) <= 0 for all I in set Z.
///ZZ() 		An M-array of working space.
///Index()		An integer working array of length at least N. 
///				On exit the contents of this array define the sets P and Z as follows...
///				Index(1) 	through	Index(NSETOP)	= Set P.
///				Index(IZ1)	through	Index(IZ2)		= Set Z.
///				IZ1 = NSETP + 1 = NPP1	
///				IZ2 = N
///Mode 		This is a success-failure flag with the following meanings.
///				1	The solution has been computed successfully.
///				2	The dimensions of the problem are bad. Either M <= 0 or N <= 0.
///				3	Iteration count exceeded. More than 3*N iterations.		

int GetIndex(int aRow, int aCol, int aNumberOfColumn)
{
	int Index;
	Index = (aRow - 1) * aNumberOfColumn + aCol;
	return Index;
}

//The following block of code is used as an internal subroutine to solve the triangular system.
//Putting the solution in ZZ[].
void Triangular(int NSETP, double ZZ[], int JJ, int Index[], double A[], int aNumberOfColumn)
{
	int II, IP, L;
	for(L = 1; L <= NSETP; L++)
	{
		IP = NSETP + 1 - L;
		if(L != 1)
		{
			for(II = 1; II <= IP; II++)
				ZZ[II] = ZZ[II] - A[GetIndex(JJ, II, aNumberOfColumn)] * ZZ[IP + 1];
		}
		JJ = Index[IP];
		ZZ[IP] = ZZ[IP] / A[GetIndex(JJ, IP, aNumberOfColumn)];
	}
}
		
void nnls(double* A, int MDA, int M, int N, double* B, double* X, double* RNorm, double* W, double* ZZ, int* Index, int* Mode)
{
	const int Zero = 0;
	const int Two = 2;
	const double Factor = 0.01f;
	int Iter, ItMax, I, J, L, IZ, IZ1, IZ2, IZMax, JZ, NSETP, NPP1, IP, II, JJ;
	double SM, WMax, Max, ASave, Up, ZTest, Alpha, CC, SS, UNorm, T;
	double* Dummy;
	Boolean Flag, TestDual, TestX, MainLoop, SecondaryLoop;
		
	// Parameter adjustments 
	--A;	
	--B;
	--X;
	--W;
	--ZZ;
	--Index;
	
	//Initialize
	Dummy = (double *)malloc(sizeof(double) * M);	
	MainLoop = TRUE;
	JJ = IZMax = 0;	
	
	(*Mode) = 1;
	if ((M <= 0) || (N <= 0))
	{
		(*Mode) = 2;
		return;
	}
	
	Iter = 0;
	ItMax = 10 * N;


	//Normalize A and B,. Divide every elemtent in A and B by maximum(all the elements in A and B) to avoid overflow problem.
	//Find Max
	Max = A[1];
	for (I = 2; I <= M * N; I++)
	{
		if (A[I] > Max)
			Max = A[I];
	}
	for (I = 1; I <= M; I++)
	{
		if (B[I] > Max)
			Max = B[I];
	}
	//Normalize
	for (I = 1; I <= M * N; I++)
		A[I] = A[I] / Max;
	for (I = 1; I <= M; I++)
		B[I] = B[I] / Max;
	
	//Initialize X
	for (I = 1; I <= N; I++)
	{
		X[I] = 0.0;
		Index[I] = I;
	}
	
	IZ2 = N;
	IZ1 = 1;
	NSETP = 0;
	NPP1 = 1;
	
	//******Main loop begins here******
	//Quit if all coefficients are already in the solution or if M columns of A have been triangularized.
	while (MainLoop == TRUE)
	{
		TestDual = TRUE;
		SecondaryLoop = TRUE;	
		
		if ((IZ1 > IZ2) || (NSETP >= M))
		{
			MainLoop = FALSE;
			break;
		}
		
		//Compute components of the dual (negative gradient) vector W().
		for (IZ = IZ1; IZ <= IZ2; IZ++)
		{
			J = Index[IZ];
			SM = (double)Zero;
			for (L = NPP1; L <= M; L++)
				SM = SM + A[GetIndex(J, L, M)] * B[L];
			W[J] = SM;
		}
			
		while(TestDual == TRUE)
		{
			TestDual = FALSE;
			//Find largest positive W[J]
			WMax = (double)Zero;
			for (IZ = IZ1; IZ <= IZ2; IZ++)
			{
				J = Index[IZ];
				if (W[J] > WMax)
				{
					WMax = W[J];
					IZMax = IZ;
				}
			}
		
			//If WMax <= 0, go to termination.
			//This indicates satisfaction of the kuhn-tucker conditions.
			if (WMax <= 0)
			{
				MainLoop = FALSE;
				break;
			}
				
			IZ = IZMax;
			J = Index[IZ];
		
			//The sign of W[J] is ok for J to be moved to set P.
			//Begin the transformation and check new diagonal element to avoid near linear dependence.
			ASave = A[GetIndex(J, NPP1, M)];
					
			H12(1, NPP1, NPP1 + 1, M, A + GetIndex(J, 1, M) - 1, 1, &Up, Dummy, 1, 1, 0);
					
			UNorm = (double)Zero;
			if (NSETP != 0)
			{
				for(L = 1; L <= NSETP; L++)
					UNorm = UNorm + (double)pow(A[GetIndex(J, L, M)], 2);
			}
			UNorm = (double)sqrt(UNorm);
		
			if (dDiff(UNorm + dAbs(A[GetIndex(J, NPP1, M)]) * Factor, UNorm) <= 0)
				Flag = TRUE;
			else
			{
				
				Flag = FALSE;
				for(L = 1; L <= M; L++)
					ZZ[L] = B[L];
								
				H12(2, NPP1, NPP1 + 1, M, A + GetIndex(J, 1, M) - 1, 1, &Up, ZZ, 1, 1, 1);
				ZTest = ZZ[NPP1] / A[GetIndex(J, NPP1, M)];
							
				//See if ZTest is positive
				if (ZTest <= 0)
					Flag = TRUE;
			}
		
			//Reject J as a candidate to be moved from set Z to set P.
			//Restore A[NPP1, J]. Set W[J] = 0. And Loop back to test dual coeffs again.
			if (Flag == TRUE)
			{
				A[GetIndex(J, NPP1, M)] = ASave;
				W[J] = 0;
				TestDual = TRUE;
			}
		}
		
		if(MainLoop != FALSE)
		{
			//The Index J = Index[IZ] has been selected to be moved from set Z to set P. Update B. Update Indices.
			//Apply householder transformations to cols in new set Z. Zero subdiagonal elts in col J. Set W[J] = 0.
			for(L = 1; L <= M; L++)
				B[L] = ZZ[L];
			
			Index[IZ] = Index[IZ1];
			Index[IZ1] = J;
			IZ1++;
			NSETP = NPP1;
			NPP1++;
		
			if(IZ1 <= IZ2)
			{
				for(JZ = IZ1; JZ <= IZ2; JZ++)
				{
					JJ = Index[JZ];
					H12(2, NSETP, NPP1, M, A + GetIndex(J, 1, M) - 1, 1, &Up, A + GetIndex(JJ, 1, M) - 1, 1, MDA, 1);
				}
			}
		
			if(NSETP != M)
			{
				for(L = NPP1; L <= M; L++)
					A[GetIndex(J, L, M)] = (double)Zero;
			}
		
			W[J] = (double)Zero;
		
			//Solve the triangular system.
			//Store the solution temporarily in ZZ[].
			Triangular(NSETP, ZZ, JJ, Index, A, M); 
		
			//******Secondary loop begins here******
			//Iteration counter.
			while(SecondaryLoop == TRUE)
			{
				TestX = TRUE;
				Iter++;			
				
				if(Iter > ItMax)
				{
					(*Mode) = 3;
					printf("Iteration exceed.");
					MainLoop = FALSE;
					break;
				}
			
				//See if all new constrained coeffs are feasible. If not compute Alpha.
				Alpha = (double)Two;
				for(IP = 1; IP <= NSETP; IP++)
				{
					L = Index[IP];
					if(ZZ[IP] <= 0)
					{					
						T = -X[L] / (ZZ[IP] - X[L]);
						if(Alpha > T)
						{
							Alpha = T;
							JJ = IP;
						}
					}
				}
		
				//If all new constrained coeffs are feasible then Alpha will still = 2. 
				//If so exit from secondary loop to main loop.
				//Otherwise use Alpha which will be between 0 and 1 to interpolate between the old X and the new ZZ.
				if(Alpha != Two)
				{
					for(IP = 1; IP <= NSETP; IP++)
					{
						L = Index[IP];
						X[L] = X[L] + Alpha * (ZZ[IP] - X[L]);
					}
					
					//Modify A and B and the index arrays to move coefficient I from set P to set Z.
					I = Index[JJ];
					while(TestX == TRUE)
					{
						X[I] = (double)Zero;
						if(JJ != NSETP)
						{
							JJ++;
							for(J = JJ; J <= NSETP; J++)
							{
								II = Index[J];
								Index[J - 1] = II;
								G1(&A[GetIndex(II, J - 1, M)], &A[GetIndex(II, J, M)], &CC, &SS, &(A[GetIndex(II, J - 1, M)]));
								A[GetIndex(II, J, M)] = (double)Zero;
								for(L = 1; L <= N; L++)
									if(L != II)
										G2(&CC, &SS, &(A[GetIndex(L, J - 1, M)]), &(A[GetIndex(L, J, M)]));
										
								G2(&CC, &SS, &(B[J-1]), &(B[J]));
							}			
						}
						
						NPP1 = NSETP; 
						NSETP--;
						IZ1--;
						//printf("IZ1 = %d", IZ1);
						Index[IZ1] = I;
						
						TestX = FALSE;
						//See if the remaining coeffs in set P are feasible.
						//They should be because of the way Alpha was determined.
						//If any are infeasible it is due to round-off error.
						//Any that are nonpositive will be set to zero and moved from set P to set Z.						
						for(JJ = 1; JJ <= NSETP; JJ++)
						{
							I = Index[JJ];
							if(X[I] < 0) {
								TestX = TRUE;
								break;
							}
						}										
					}
					
					//Copy B[] into ZZ[] 
					//Then solve again and loop back.
					for(I = 1; I <= M; I++)
						ZZ[I] = B[I];
					
					Triangular(NSETP, ZZ, JJ, Index, A, M); 
				}
				else
					SecondaryLoop = FALSE;			
			}
			
			//******End of secondary loop******
			if(SecondaryLoop == FALSE)
			{
				for(IP = 1; IP <= NSETP; IP++)
				{
					I = Index[IP];
					X[I] = ZZ[IP];
				}	
			}
			//All new coeffs are positive. Loop back to main loop.		
			
		}
	}

	//******End of Main loop******
	//Come to here to termination.
	//Compute the norm of the final residual vector.
	SM = (double)Zero;
	if (NPP1 > M)
	{
		for (J = 1; J <= N; J++)
			W[J] = (double)Zero;
			
		(*RNorm) = (double)sqrt(SM);
	}
	else
	{
		for (I = NPP1; I <= M; I++)
			SM = SM + (double)pow(B[I], 2);
			
		(*RNorm) = (double)sqrt(SM);
	}
	return;	
}

#include <time.h>
#include "nnls.h"

int main()
{
	int MDA;
	int Mode;
	double RNorm;
	const int M = 100;
	const int N = 100;
	int i;	
	double *A, *W, *ZZ, *B, *X;
	int *Index;
	FILE *fp;

	printf("test.c: In test Main Start.\n");
	
    //Test double version Start
	printf("test.c: Test double Start.\n");
	fflush(stdout);

	Mode = 0;
	RNorm = 0;
	MDA = M;
	A = (double *)malloc(sizeof(double) * M * N);
	W = (double *)malloc(sizeof(double) * N);
	ZZ = (double *)malloc(sizeof(double) * M);
	B = (double *)malloc(sizeof(double) * M);
	X = (double *)malloc(sizeof(double) * N);
	Index = (int *)malloc(sizeof(int) * N);
				
	if (!A || !W || !ZZ || !B || !X || !Index)
	{
		printf("test.c: Malloc error.\n");
		return 1;
	}
		
	printf("test.c: Initialize.\n");
	fflush(stdout);
	
	srand( (unsigned)time( NULL ) );	
	for(i = 0; i < (M * N); i++)
	{
		A[i] = rand() % 50;
	}
		
	if((fp=fopen("C.txt", "w+"))==NULL) 
	{
		printf("test.c: Cannot open C.txt.\n");
		exit(1);
	}

	for(i = 0; i < (M * N); i++)
	{
		fprintf(fp, "%lf", A[i]);
		if (i % 100 == 99)
			fprintf(fp, "\n");
		else
			fprintf(fp, " ");
	}
	fclose(fp);	
	
	srand( (unsigned)time( NULL ) );
	for(i = 0; i < M; i++)
	{

		B[i] = rand() % 100;
	}
	
	if((fp=fopen("D.txt", "w+"))==NULL) 
	{
		printf("test.c: Cannot open D.txt.\n");
		exit(1);
	}

	for(i = 0; i < M; i++)
	{
		fprintf(fp, "%lf ", B[i]);
	}
	fclose(fp);
	
	for(i = 0; i < N; i++)
		X[i] = 0;
	fflush(stdout);	
	nnls(A, MDA, M, N, B, X, &RNorm, W, ZZ, Index, &Mode);

	if((fp=fopen("Xf.txt", "w+"))==NULL) 
	{
		printf("test.c: Cannot open Xf.txt.\n");
		exit(1);
	}

	for(i = 0; i < N; i++)
	{
		fprintf(fp, "%8.7e ", X[i]);
	}
	fclose(fp);	

	free(A);
	free(W);
	free(ZZ);
	free(B);
	free(X);
	free(Index);
	printf("test.c: Test NNLS End.\n");
//Test double end
	
	printf("test.c: Main end.\n");
	return 0;
}

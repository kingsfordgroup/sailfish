/*
 * MLR-MCL (Multi-Level Regularized Markov Clustering) - Version 1.2
 * 
 * Copyright 2010, The Ohio State University. 
 * Portions Copyright 1997, Regents of the University of Minnesota.
 * All rights reserved.
 * 
 * Redistribution and use in source and binary forms, with or
 * without modification, are
 * permitted provided that the following conditions are met:
 * 
 * 1. Redistributions, with or without modifications, of source 
 * code must retain the above
 * copyright notice, this list of conditions and the following
 * disclaimer.
 * 
 * 2. Redistributions, with or without modifications, in binary form 
 * must reproduce the above
 * copyright notice, this list of conditions and the following
 * disclaimer
 * in the documentation and/or other materials provided with the
 * distribution.
 * 
 * 3. The names of the Ohio
 * State University, the University of Minnesota and
 * their contributors may not be used to endorse or promote products
 * derived from this software without specific prior permission. 
 * 
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND
 * CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
 * FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 * COPYRIGHT
 * OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
 * INCIDENTAL,
 * SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
 * LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
 * USE,
 * DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON
 * ANY
 * THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
 * TORT
 * (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE
 * USE
 * OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH
 * DAMAGE.
 * 
 */

/*
 * Contributors: Venu Satuluri, Srinivasan Parthasarathy
 * 
 * Reference: "Scalable Graph Clustering using Stochastic Flows",
 * KDD 2009. http://doi.acm.org/10.1145/1557019.1557101
 */


/* This file has the core mcl routines
 */

/* 
 * $Id $
 */
 #include <metis.h>

void swapwgttype(wgttype* a, wgttype* b)
{
	wgttype t;
	t=*a;
	*a=*b;
	*b=t;
}

void swapidxtype(idxtype* a, idxtype* b)
{
	idxtype t;
	t=*a;
	*a=*b;
	*b=t;
}

int insertLastElement(idxtype* adjncy, wgttype* adjwgt, int count)
{
	/* This function assumes adjncy is an array of size count 
	   that is sorted in
	 * increasing order, except for the last element. 
	 * this function will do a series of swaps beginning from
	 * the last but one element until it is in the right place in
	 * the array. adjwgt will have the same swaps performed on it
	 * as those on adjncy.
	 */	 
	 int i;
	 for( i=count-1; i>0 && (adjncy[i] < adjncy[i-1]) ; i--)
	 {
		idxtype t1;
	 	wgttype t2;

		t1=adjncy[i];
		adjncy[i]=adjncy[i-1];
		adjncy[i-1]=t1;

		t2=adjwgt[i];
		adjwgt[i]=adjwgt[i-1];
		adjwgt[i-1]=t2;
	 	
//		swapidxtype((adjncy+i),(adjncy+i-1));
//	 	swapwgttype((adjwgt+i),(adjwgt+i-1));
	 }
	 return i;
}

int compareints(const void * a, const void * b)
{
	return ( *(int*)a - *(int*)b );
}

int comparewgttypes(const void * a, const void * b)
{
	// need the array in descending order!
	wgttype diff= (*(wgttype*)b - *(wgttype*)a );
	return ( (diff>0) ? 1 : ((diff<0) ? -1: 0) );
}
  
wgttype computeThreshold(wgttype avg, wgttype max)
{
	wgttype ret = MLMCL_PRUNE_A*avg*(1-MLMCL_PRUNE_B*(max-avg));	
//	wgttype ret = 0.8*avg*(1-4*(max-avg));	
//	wgttype ret = (wgttype)pow((double)max,3.0);
//	return (ret>0.00001)?ret:0.00001; 
	ret = (ret>1.0e-7)?ret:1.0e-7; 
	ret = (ret > max) ? max : ret;

	return ret; 
//	return 0;
}

wgttype threshPrune(int *n, idxtype* adjncy, wgttype* adjwgt,
wgttype thresh, wgttype initsum)
{
	idxtype*
	indicesToRetain=(idxtype*)malloc(sizeof(idxtype)*(*n));
	int i,j;
	wgttype sum=0;
	
	for(i=0,j=0;i<*n;i++)
	{
		if ( adjwgt[i]>=thresh )
		{
			sum+=adjwgt[i];
			indicesToRetain[j++]=i;
		}
	}
	if ( j == 0 )
	{
/*		free(indicesToRetain);
		printf("Yikes! new_vector_nnz=0!\n");
		return initsum; */
//		printf("Threshold:%f\n", thresh);
//		abort();
	}
	for(i=0;i<j;i++)
	{
		adjncy[i]=adjncy[indicesToRetain[i]];
		adjwgt[i]=adjwgt[indicesToRetain[i]];
	}
	*n=j;
	free(indicesToRetain);
	return sum;
	
}

void exactPruneGraph(GraphType* graph, int l, idxtype thresh )
{
/*	Matrix M;
	int i;
	M.nnz = graph->nedges;
	M.adjncy = graph->adjncy;
	M.adjwgt = (wgttype*) malloc( sizeof(wgttype)*M.nnz );
	for ( i=0; i<M.nnz; i++)
		M.adjwgt[i] = (wgttype)graph->adjwgt[i];

	exactPruneMatrix(&M, l);

	graph->nedges = M.nnz;
	graph->adjncy = M.adjncy;
	for ( i=0; i<M.nnz; i++)
		graph->adjwgt[i] = (idxtype) M.adjwgt[i];
	graph->adjwgt = (idxtype*) realloc (graph->adjwgt, sizeof(wgttype)*M.nnz);

	return;
*/
	int i, j,k;
	long nnz = graph->nedges;

	if ( thresh == 0 )
	{
		long s = sizeof(idxtype)*(nnz);
		idxtype* newadjwgt=(idxtype*)malloc(s);

		for(i = 0; i < nnz; i++)
			newadjwgt[i]=graph->adjwgt[i];
			
		thresh = RandomSelectInts(newadjwgt, 0, nnz-1, nnz-l+1);
		free(newadjwgt);
	}

//	printf("Threshold:%d\n",thresh);

	for(i=0,j=0; i<graph->nvtxs; i++)
	{
		k=graph->xadj[i];
//		if ( i > 0 )
//		{
			graph->xadj[i]=j;
//			printf("M->xadj[%d] set to %d\n", i, j);
//			printf("M->xadj[%d] is %d\n",i+1,M->xadj[i+1]);
//		}

		for(; k<graph->xadj[i+1]; k++)
		{
//			printf("M->adjwgt[%d]:%.2f\n",k,M->adjwgt[k]);
			if ( graph->adjwgt[k] >= thresh 
				|| graph->adjncy[k]	== i )
			{
				if ( j < k )
				{
					graph->adjncy[j]=graph->adjncy[k];
					graph->adjwgt[j]=graph->adjwgt[k];
				}
				j++;
			}
		}
//		printf("j:%d\n",j);
	}
	graph->xadj[graph->nvtxs]=j;
//	printf("M->nnz set to %d\n", j);
	
	if ( j < graph->nedges )
	{
		graph->nedges = j;
		graph->adjncy = (idxtype*) realloc(graph->adjncy,j*sizeof(idxtype));
		graph->adjwgt =	(idxtype*) realloc(graph->adjwgt,j*sizeof(idxtype));
	}

}

/* Retains the top l entries in the sparse matrix M */
void exactPruneMatrix(Matrix* M, int l)
{
	wgttype thresh;
	int i,j,k;
	wgttype* newadjwgt=(wgttype*)malloc(sizeof(wgttype)*M->nnz);

	for(i=0;i<M->nnz;i++)
		newadjwgt[i]=M->adjwgt[i];
		
//	printf("Pruning matrix with l:%d\n", l);

	thresh=RandomSelect(newadjwgt, 0, M->nnz-1,	M->nnz-l+1);
	free(newadjwgt);
//	dumpMatrix(M);

//	printf("Threshold:%f\n",thresh);

	for(i=0,j=0; i<M->nvtxs; i++)
	{
		k=M->xadj[i];
//		if ( i > 0 )
//		{
			M->xadj[i]=j;
//			printf("M->xadj[%d] set to %d\n", i, j);
//			printf("M->xadj[%d] is %d\n",i+1,M->xadj[i+1]);
//		}

		for(; k<M->xadj[i+1]; k++)
		{
//			printf("M->adjwgt[%d]:%.2f\n",k,M->adjwgt[k]);
			if ( M->adjwgt[k] >= thresh )
			{
				if ( j < k )
				{
					M->adjncy[j]=M->adjncy[k];
					M->adjwgt[j]=M->adjwgt[k];
				}
				j++;
			}
		}
//		printf("j:%d\n",j);
	}
	M->xadj[M->nvtxs]=j;
//	printf("M->nnz set to %d\n", j);
	
	if ( j < M->nnz )
	{
		M->nnz=j;
		M->adjncy=(idxtype*)realloc(M->adjncy,j*sizeof(idxtype));
		M->adjwgt=(wgttype*)realloc(M->adjwgt,j*sizeof(wgttype));
	}

}


/*
 * Does exact pruning of one column. n indicates number of column
 * entries, adjncy is an array of node ids, adjwgt is the array
 * of flow values and k is the number of entries that need to be
 * retained. The return value is the sum of the retained entries.
 */
wgttype exactPrune(int *n, idxtype* adjncy, wgttype* adjwgt,
						int k, wgttype initsum)
{
	int i,j,initn=*n;
	wgttype sum=0, *newadjwgt,thresh;
	// if column size less than k, don't need to do anything.
	if ( k >= *n )
		return initsum;

	// RandomSelect moves elements in the array around. Hence,
	// need a placeholder array.
	newadjwgt=(wgttype*)malloc(sizeof(wgttype)*(*n));
	for( i=0; i<*n; i++ )
	{
		newadjwgt[i]=adjwgt[i];
	}

/*	qsort((void*)newadjwgt, *n, sizeof(wgttype),
	comparewgttypes); */
//	thresh=newadjwgt[k-1];

	thresh=RandomSelect(newadjwgt,0,*n-1,*n-k+1);
//	printf("The %d biggest element selected by exactPrune",k);
//	printf(" is %f\n", thresh);
	free(newadjwgt);

	for(i=0,j=0;i<*n;i++)
	{
		if  ( adjwgt[i] >= thresh )
		{
			if ( j < i )
			{
				adjncy[j]=adjncy[i];
				adjwgt[j]=adjwgt[i];
			}
			sum += adjwgt[i];
			j++;
		}
	}
	*n=j;

/*	if ( *n > k )
	{
		printf("*n > k! *n:%d, k:%d, thresh:%1.4f\n",*n,k,thresh);
		printf("%1.1f",newadjwgt[0]);
		for(i=1;i<*n;i++)
			printf(",%1.1f",newadjwgt[i]);
		printf("\n");
		abort();
	}
*/
	if (!( sum > 0)  || sum > 1 )
	{
		printf("sum:%f,n:%d,k:%d,initsum:%f,thresh:%f\n",
					sum,*n,k,initsum,thresh); 
/*		for(i=0;i<initn;i++)
			printf("%f ", newadjwgt[i]);
		printf("\n");
*/		abort();
	}

//	free(newadjwgt);
	return sum;
}

wgttype* getRowSums(Matrix* a)
{
	wgttype* rowSums = (wgttype*) malloc(sizeof(wgttype)*a->nvtxs); 

	int i;

	for ( i=0; i<a->nvtxs; i++ )
		rowSums[i] = 0;

	for ( i=0; i<a->xadj[a->nvtxs]; i++ )
		rowSums[a->adjncy[i]]+=a->adjwgt[i];

	return rowSums;
}

Matrix* getTranspose2(Matrix* M)
{
	int nvtxs = M->nvtxs;
	Matrix* ret=allocMatrix(M->nvtxs, M->xadj[nvtxs], 0, 0, 0);
	idxtype *inDegrees = idxmalloc(M->nvtxs,
	"getTranspose2:inDegrees");
	for ( int i=0; i<nvtxs; i++ )
		inDegrees[i] = 0;
	for ( int i=0; i<M->xadj[nvtxs]; i++ )
		inDegrees[M->adjncy[i]]++;

	ret->xadj[0]=0;
	for ( int i = 0; i < nvtxs; i++)
		ret->xadj[i+1] = ret->xadj[i] + inDegrees[i];
	
	idxtype *counters = inDegrees;
	for ( int i = 0; i < nvtxs; i++)
		counters[i] = ret->xadj[i];

	for ( int i = 0; i < M->nvtxs; i++ )
	{
		for ( int j = M->xadj[i]; j < M->xadj[i+1]; j++ )
		{
			int k = M->adjncy[j];
			ret->adjncy[counters[k]] = i;
			ret->adjwgt[counters[k]] = M->adjwgt[j];
			counters[k]++;
		}
	}
	return ret;
}

Matrix* getTranspose(Matrix* M)
{
	return getTranspose2(M);
	//printf("M->nvtxs:%d, M->nnz:%d\n",M->nvtxs,M->nnz);
	Matrix* ret=allocMatrix(M->nvtxs, M->nnz, 0, 0, 0);
	idxtype **temp_adjncy, *adjncy=M->adjncy, nvtxs=M->nvtxs,
	nedges=M->nnz, *xadj=M->xadj;
	wgttype **temp_adjwgt, *adjwgt=M->adjwgt;
	int* currentSizes=(int*)malloc(sizeof(int)*nvtxs);
	int* allocSizes=(int*)malloc(sizeof(int)*nvtxs);
	int i,j,initSize=nedges/nvtxs;
	initSize=(initSize>20)?initSize:20;

	temp_adjncy=(idxtype**)malloc(sizeof(idxtype*)*nvtxs);
	temp_adjwgt=(wgttype**)malloc(sizeof(wgttype*)*nvtxs);

	for( i=0; i<nvtxs; i++)
	{
		currentSizes[i]=0;
		allocSizes[i]=initSize;
		temp_adjncy[i]=idxmalloc(initSize, "temp_adjncy");
		temp_adjwgt[i]=(wgttype*)malloc(initSize*sizeof(wgttype));
	}

	for(i=0;i<nvtxs;i++)
	{
		for(j=xadj[i];j<xadj[i+1];j++)
		{
			idxtype from=i,to=adjncy[j];
			wgttype wgt;
			wgt=adjwgt[j];
			currentSizes[to]++;
			if ( currentSizes[to] > allocSizes[to] )
			{
				allocSizes[to]+=allocSizes[to];
				if ( allocSizes[to] == 0 )
				{
					printf("allocSizes[%d] is 0, initSize:%d\n",
								to, initSize);
					abort();
				}
				temp_adjncy[to]=(idxtype*)realloc(temp_adjncy[to],allocSizes[to]*sizeof(idxtype));
				if ( temp_adjncy[to] == NULL )
				{
					printf("Could not allocate %ld",
							allocSizes[to]*sizeof(idxtype)); 
					printf(" bytes in getTranspose");
					abort();
				}
				temp_adjwgt[to]=(wgttype*)realloc(temp_adjwgt[to],allocSizes[to]*sizeof(idxtype));
				if ( temp_adjwgt[to] == NULL )
				{
						printf("Could not allocate %ld",
								allocSizes[to]*sizeof(idxtype)); 
						printf(" bytes in getTranspose");
						abort();
				}
			}
			temp_adjncy[to][currentSizes[to]-1]=from;
			temp_adjwgt[to][currentSizes[to]-1]=wgt;
		}
	}

	ret->xadj[0]=0;
	for( i=0; i<nvtxs; i++ )
	{
		ParallelQSort(temp_adjncy[i],
				temp_adjwgt[i],0,currentSizes[i]-1);

		ret->xadj[i+1]=ret->xadj[i]+currentSizes[i];
		for(j=ret->xadj[i]; j<ret->xadj[i+1]; j++)
		{
			ret->adjncy[j]=temp_adjncy[i][j-ret->xadj[i]];
			ret->adjwgt[j]=temp_adjwgt[i][j-ret->xadj[i]];
		}

		free(temp_adjncy[i]);
		free(temp_adjwgt[i]);
	}

	free(temp_adjncy);
	free(currentSizes);
	free(allocSizes);
	free(temp_adjwgt);

	return ret;
}

Matrix* expand(Matrix* M, Matrix* M0)
{
	wgttype *M_adjwgt, *M0_adjwgt,wgt;
	idxtype *M_xadj, *M_adjncy, *M0_adjncy, *M0_xadj;
	int retM_size, i, j, k, l, *index, nnzcount=0, M_nedges;
	timer bsearchTimer;
	cleartimer(bsearchTimer);
	
	Matrix* ret=allocMatrix(M->nvtxs, M->nnz,1,1,0);
	/*allocate ret same amount of nonzeros as M. */
	ret->nvtxs=M->nvtxs;
	retM_size=M_nedges=M->nnz;

	M_xadj=M->xadj;
	M_adjncy=M->adjncy;
	M_adjwgt=M->adjwgt;
	M0_xadj=M0->xadj;
	M0_adjncy=M0->adjncy;
	M0_adjwgt=M0->adjwgt;
	
	ret->xadj[0]=0;
	for( i=0; i<M->nnz; i++)
		ret->adjncy[i]=-1;

/*	for( i=0; i<M->nnz;i++)
		printf("M_adjwgt[i]:%1.1f\n",M_adjwgt[i]);
*/
/*	for( i=0; i<M->nvtxs;i++ )
		printf("M->adjwgtsum[%d]:%1.4f\n",i,M->adjwgtsum[i]);
*/
//	dumpMatrix(M);
//	dumpMatrix(M0);
	for( i=0; i<M->nvtxs; i++ )
	{
	 	/* iterate through each of i's neighbours from the
		 * initial adjacency matrix, i.e. M0_adjncy*/
		idxtype ibase=ret->xadj[i];
	 	for( j=M0_xadj[i]; j < M0_xadj[i+1]; j++)
		{
			/* k is the neighbour's id and wgt is its weight, that we'll use to
			 * weight its column */
			k=M0_adjncy[j];
			wgt=M0_adjwgt[j];
			for( l=M_xadj[k]; l<M_xadj[k+1]; l++ )
			{
				starttimer(bsearchTimer);
				index =	(int*)bsearch((M_adjncy+l),
							(ret->adjncy+ibase),(nnzcount-ibase),
							sizeof(int), compareints);
				stoptimer(bsearchTimer);
				if ( index == NULL )
				{
					idxtype kk;
					nnzcount++;
					if ( nnzcount > retM_size )
					{
						ret->adjncy=(idxtype*)realloc(ret->adjncy,
						(retM_size+M_nedges)*sizeof(idxtype));
						ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
						(retM_size+M_nedges)*sizeof(wgttype));
						retM_size+=M_nedges;
					}	
					ret->adjncy[nnzcount-1]=M_adjncy[l];
					ret->adjwgt[nnzcount-1]=(wgttype)wgt*M_adjwgt[l];
			/*		if ( ret->adjwgt[nnzcount-1] == 0 )
					{
						printf("before insertLastElement,");
						printf("adjwgt:%1.2f, wgt:%1.2f,",
						ret->adjwgt[nnzcount-1],wgt);
						printf("M_adjwgt[j]:%1.2f\n",M_adjwgt[l]);
						abort();
					}
			*/		kk=insertLastElement((ret->adjncy+ibase),
								(ret->adjwgt+ibase),nnzcount-ibase);
					
/*					printf("%d inserted into list of %d	at position %d",
							M_adjncy[l], i, kk);
					printf(" with wgt:%1.4f\n",ret->adjwgt[nnzcount-1]);
					printf("nbr wgt:%1.4f, nbr flow:%1.4f\n",
							wgt, M_adjwgt[l]);
*/
				/*	if ( ret->adjwgt[nnzcount-1] == 0 )
					{
						printf("adjwgt 0,");
						printf("wgt:%1.2f,", wgt);
						printf("M_adjwgt[j]:%1.2f\n",M_adjwgt[l]);
						abort();
					} */
				}
				else
				{
					//nnzcount++;
					wgttype	oldwgt=ret->adjwgt[ibase+(index-ret->adjncy-ibase)];
					wgttype newwgt;

/*					printf("index for %d in %d is %d\n",M_adjncy[l],
							i,(index-ret->adjncy-ibase));
*/					newwgt=oldwgt+wgt*M_adjwgt[l];
					ret->adjwgt[ibase+(index-ret->adjncy-ibase)]=newwgt;

/*					printf("wgt updated from %1.4f to %1.4f\n",
							oldwgt, newwgt);
					printf("nbr wgt:%1.4f, nbr flow:%1.4f\n",
							wgt, M_adjwgt[l]);
*/				}
				
			}
		}
		ret->xadj[i+1]=nnzcount;
		ret->nnz=nnzcount;
	 }

	// dumpMatrix(ret);
	/* for( i=0; i<ret->nvtxs;i++)
		printf("ret_adjwgtsum[i]:%1.5f\n",ret->adjwgtsum[i]);
	 */
	 printf("Time for bsearch:%7.3f\n", gettimer(bsearchTimer));
	 return ret;
}

/*
 * This function computes
 * outerDiagonalFactor*A^{T}*innerDiagonalFactor*A*outerDiagonalFactor
 * outerDiagonalFactor and innerDiagonalFactor are diagonal
 * matrices, and hence only the diagonal entries are given in a
 * array each. If these arrays are null, A^{T}A is computed.
 *
 * The resulting matrix is going to be computed block-wise, with
 * the size of the block given in blockSize. If
 * A.nvtxs=k*blockSize, then the resulting matrix is going to be
 * computed in k^2 blocks.
 * Let
 */
Matrix* AtimesATransposeBlocking(Matrix* A, wgttype* outerDiagonalFactor,
	wgttype* innerDiagonalFactor, int blockSize, wgttype
	threshold )
{
	
	/* sort adjacency lists in A. */
	int i;
	for ( i=0; i<A->nvtxs; i++ )
	{
		ParallelQSort(A->adjncy, A->adjwgt, A->xadj[i], A->xadj[i+1]-1);
	}

	int numColumnBlocks = (int) ceil( (float)A->nvtxs/(float)blockSize	);
	long ret_size = A->nvtxs*100;
	Matrix *ret = allocMatrix(A->nvtxs, ret_size, 0, 0, 0);
	long nnz = 0;
	ret->xadj[0] = 0;
	int numBlocksComputed = 0;

	timer sortTimer, totalTimer;
	cleartimer(sortTimer);
	cleartimer(totalTimer);

	starttimer(totalTimer);

	for ( i=0; i<numColumnBlocks; i++ )
	{
		int j, rowStart, rowEnd, colStart, colEnd;
		colStart = i*blockSize;
		colEnd = colStart + blockSize - 1;
		colEnd = (colEnd > A->nvtxs-1) ? A->nvtxs-1 : colEnd;

		int numRowBlocks = numColumnBlocks - i;
		Matrix** blocks = (Matrix**) 
					malloc(sizeof(Matrix*) * numRowBlocks);
		long colBlock_nnz = 0;
		for ( j=0; j<numRowBlocks; j++ )
		{
			rowStart = colStart + j*blockSize;
			rowEnd = rowStart + blockSize - 1;
			rowEnd = (rowEnd > A->nvtxs-1) ? A->nvtxs-1 : rowEnd;
			blocks[j] = AtimesATranspose(A, outerDiagonalFactor,
				innerDiagonalFactor, rowStart, rowEnd, colStart,
				colEnd, threshold, &sortTimer);
			colBlock_nnz += blocks[j]->nnz;
			numBlocksComputed++;
			if ( numBlocksComputed % 50 == 0 )
			{
				//printf("Computed %d blocks,", numBlocksComputed);
				//printf("nnz:%ld, sorting time:%.3f\n", nnz+colBlock_nnz,
				//	gettimer(sortTimer) );
				fflush(stdout);
			}
		}

		if ( nnz + colBlock_nnz > ret_size )
		{
			ret_size += colBlock_nnz;
			ret->adjncy = (idxtype*) 
					realloc( ret->adjncy, sizeof(idxtype)* ret_size );
			ret->adjwgt = (wgttype*)
					realloc( ret->adjwgt, sizeof(wgttype)* ret_size );
		}

		int k,l;
		for ( k=colStart; k <= colEnd; k++ )
		{
			for ( j=0; j<numRowBlocks; j++ )
			{
				rowStart = colStart + j*blockSize;
				rowEnd = rowStart + blockSize - 1;
				rowEnd = (rowEnd > A->nvtxs-1) ? A->nvtxs-1 : rowEnd;

				for ( l = blocks[j]->xadj[k-colStart];
						l < blocks[j]->xadj[k-colStart+1]; l++ )
				{
					int to = blocks[j]->adjncy[l] + rowStart;
					if ( to >= A->nvtxs )
					{
						printf("to:%d, j:%d, l:%d, rowStart:%d\n"
								,to, j, l, rowStart);
						abort();
					}
					ret->adjncy[nnz] = to;
					ret->adjwgt[nnz] = blocks[j]->adjwgt[l];
					nnz++;
				}
			}
			ret->xadj[k+1] = nnz;
		}

		for ( j=0; j<numRowBlocks; j++ )
		{
			freeMatrix( blocks[j] );
		}
	}
	ret->nnz = nnz;

	Matrix *retT = getTranspose(ret);
	Matrix *result = add(ret, retT); 

	stoptimer( totalTimer );
	//printf("Total time on blocking matrix multiply: %7.3f\n", 
	//			gettimer(totalTimer) );
	//printf("Time spent sorting:%7.3f\n", gettimer(sortTimer) );

	freeMatrix(ret);
	freeMatrix(retT);

	return result;
}

/* Computes one block of the result of AtimesATranspose. The
 * corners of the block are specified by rowStart, rowEnd,
 * colStart and colEnd (all inclusive). 
 * This is called as a subroutine by AtimesATransposeBlocking. 
 * NOTE! The adjacency lists of A must be sorted before this
 * function is called.
 */
Matrix* AtimesATranspose(Matrix* A, wgttype*
	outerDiagonalFactor, wgttype* innerDiagonalFactor, int rowStart,
	int rowEnd, int colStart, int colEnd, wgttype threshold,
	timer* sortTimer)
{
	int blockSize = colEnd-colStart+1;

	/* We'll only do something if the required block overlaps
	 * with the lower triangular portion of AA^{T}. We don't want
	 * to compute the upper triangular portion since the matrix
	 * is symmetric. */
	if ( rowEnd <= colStart )
		return NULL;

	/* The partial sums are stored as <key,value> pairs, where
	 * key is the <j,i> index of the matrix entry. <j,i> is
	 * stored as i*blockSize + j. The reversed order of i and j
	 * is intentional; it is because we'll be calculating and
	 * storing the result matrix in column-major order.*/
	int partialSums_size = blockSize*100;
	int partialSums_size_inc = blockSize*100;
	long* partialSumKeys = (long*) 
				malloc( sizeof(long) * partialSums_size);
	wgttype* partialSums = (wgttype*) 
				malloc( sizeof(wgttype) * partialSums_size );

	int i, j, k, l;
	long partialSum_count = 0;
	for ( i=0; i<A->nvtxs; i++ )
	{
		/* Only compute lower triangular part of resulting
		 * matrix. i.e. only those <k,j> partial 
	 	 * sums where j<k. */

		int jStart, jEnd, kStart, kEnd;
		if ( rowStart == 0 && colStart == 0 && rowEnd == A->nvtxs
				&& colEnd == A->nvtxs )
		{
			jStart = kStart = A->xadj[i];
			jEnd = kEnd = A->xadj[i+1];
		}
		else
		{
			jStart = bsearch_insertPos( A->adjncy, colStart,
						A->xadj[i], A->xadj[i+1]-1);
			jEnd = bsearch_insertPos( A->adjncy, colEnd,
					A->xadj[i], A->xadj[i+1] - 1);
			kStart = bsearch_insertPos( A->adjncy, rowStart,
						A->xadj[i], A->xadj[i+1] - 1);
			kEnd = bsearch_insertPos( A->adjncy, rowEnd,
						A->xadj[i], A->xadj[i+1] - 1);

			if ( jEnd < A->xadj[i+1] && A->adjncy[jEnd] == colEnd )
			{
				// can't have jEnd be greater than A->xadj[i+1]
				jEnd++;
			}
			if ( kEnd < A->xadj[i+1] && A->adjncy[kEnd] == rowEnd )
			{
				// can't have kEnd be greater than A->xadj[i+1]
				kEnd++; 
			}
		}

		for ( j=jStart; j<jEnd; j++ )
		{
			if ( A->adjncy[j] == i )
				continue;
			for ( k = (kStart<=j ? (j+1) : kStart); k < kEnd; k++ )
			{
				if ( A->adjncy[k] == i || A->adjncy[k] <= A->adjncy[j] )
				{
					continue;
				}

				partialSum_count++;
				if ( partialSum_count > partialSums_size )
				{
					partialSums_size += partialSums_size_inc;
					long s = partialSums_size*sizeof(long);
					partialSumKeys = (long*) realloc( 
						partialSumKeys, s );
					if ( partialSumKeys == NULL )
					{
						fprintf(stderr, "Could not allocate ");
						fprintf(stderr, "%ld bytes for ", s);
						fprintf(stderr, "partialSumKeys!\n");
						abort();
					}

					s = partialSums_size*sizeof(wgttype);
					partialSums = (wgttype*) realloc( 
						partialSums, s );
					if ( partialSums == NULL )
					{
						fprintf(stderr, "Could not allocate ");
						fprintf(stderr, "%ld bytes for ", s);
						fprintf(stderr, "partialSums!\n");
						abort();
					}
				}
	
				long t = partialSum_count - 1;
				int offsetAdjj = A->adjncy[j] - colStart;
				int offsetAdjk = A->adjncy[k] - rowStart;
				if ( offsetAdjj < 0 || offsetAdjk < 0 )
				{
					fprintf(stderr, "Yikes, A->adjncy[j]:%d, ",
									A->adjncy[j]);
					fprintf(stderr, "A->adjncy[k]:%d", A->adjncy[k]);
					fprintf(stderr, ", rowStart:%d", rowStart);
					fprintf(stderr, ", colStart:%d", colStart);
					abort();
				}
				partialSumKeys[t] = offsetAdjj*blockSize + offsetAdjk;
				partialSums[t] = A->adjwgt[j] * A->adjwgt[k];
				if ( innerDiagonalFactor != NULL )
					partialSums[t] *= innerDiagonalFactor[i];
				if ( outerDiagonalFactor != NULL )
				{
					partialSums[t] *=
					outerDiagonalFactor[A->adjncy[j]] *
					outerDiagonalFactor[A->adjncy[k]];
				}

			}
		}
	}

//	printf("Going to sort %d partialSums\n", partialSum_count);
	starttimer(*sortTimer);
	ParallelQSortLongs(partialSumKeys, partialSums, 0,
				partialSum_count-1);	
	stoptimer(*sortTimer);
//	printf("Done sorting\n");

	int nnz=1;
	for ( i=1; i<partialSum_count; i++ )
	{
		if ( partialSumKeys[i] != partialSumKeys[i-1] )
			nnz++;
	}

	Matrix* ret = allocMatrix(blockSize, nnz, 0, 0, 0);

	/* Do reduction of partialSums to fill in matrix. */
	ret->xadj[0] = 0;
	int lasti = 0, lastj = 0, nnzcount = 0;
	for ( k=0; k<partialSum_count; k++ )
	{
		i = partialSumKeys[k] / blockSize;
		j = partialSumKeys[k] % blockSize;

		if ( k > 0 && lasti == i && lastj == j )
		{	
			ret->adjwgt[nnzcount] += partialSums[k];
		}
		else
		{
			if ( k > 0 && ret->adjwgt[nnzcount] >= threshold )
				nnzcount++;
				
			for ( ; lasti != i; lasti++ )
			{
				ret->xadj[lasti+1] = nnzcount;
			}
			ret->adjncy[nnzcount] = j;
			ret->adjwgt[nnzcount] = partialSums[k];
		}
		lasti = i; lastj = j;
	}

	if ( ret->adjwgt[nnzcount] >= threshold )
		nnzcount++;

	for ( ; i < blockSize; i++ )
		ret->xadj[i+1] = nnzcount;
	ret->nnz = nnzcount;
	if ( nnz > nnzcount )
	{
		ret->adjwgt = (wgttype*)realloc( ret->adjwgt, 
						sizeof(wgttype)*nnzcount );
		ret->adjncy = (idxtype*)realloc( ret->adjncy,
						sizeof(idxtype)*nnzcount );
	}

	free ( partialSumKeys );
	free ( partialSums );

	return ret;
}

Matrix* expand_ht(Matrix* M, Matrix* M0, idxtype*
	hashtable, int sortLists, wgttype threshold)
{
	wgttype *M_adjwgt, *M0_adjwgt,wgt;
	idxtype *M_xadj, *M_adjncy, *M0_adjncy, *M0_xadj, htcounter;
	int i, j, k, l, index, nnzcount=0, M_nedges;
	int init_retM_size;
	long retM_size;

	wgttype max = 0, sum = 0, min = 0;
	
//	printf("nvtxs:%d, nnz:%d\n",M->nvtxs, M->nnz);
	Matrix* ret=allocMatrix(M0->nvtxs, M->nnz,1,1,0);
	/*allocate ret same amount of nonzeros as M. */
	ret->nvtxs=M0->nvtxs;
	retM_size=M_nedges=M->nnz;
	init_retM_size = retM_size;

	M_xadj=M->xadj;
	M_adjncy=M->adjncy;
	M_adjwgt=M->adjwgt;
	M0_xadj=M0->xadj;
	M0_adjncy=M0->adjncy;
	M0_adjwgt=M0->adjwgt;
	
	for(i=0;i<M0->nvtxs;i++)
	{
		hashtable[i]=-1;
	}

	ret->xadj[0]=0;
	nnzcount=0;
	long numPrunedAway = 0;
	int numZeroCols = 0;
	for( i=0; i<M0->nvtxs; i++ )
	{
	 	/* iterate through each of i's neighbours from the
		 * initial adjacency matrix, i.e. M0_adjncy*/
		idxtype ibase=ret->xadj[i];
		htcounter=0;
	 	for( j=M0_xadj[i]; j < M0_xadj[i+1]; j++)
		{
			/* k is the neighbour's id and wgt is its weight, that we'll use to
			 * weight its column */
			k=M0_adjncy[j];
			wgt=M0_adjwgt[j];
/*			if ( i == 250 && k == 740 )
			{
				printf("Yay!\n");
			}
*/
			for( l=M_xadj[k]; l<M_xadj[k+1]; l++ )
			{
				index=hashtable[M_adjncy[l]];
				if ( index == -1 )
				{
					nnzcount++;
					if ( nnzcount > retM_size )
					{
						retM_size += init_retM_size;
						//retM_size+=retM_size;
						ret->adjncy=(idxtype*)realloc(ret->adjncy,
						(retM_size)*sizeof(idxtype));
						ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
						(retM_size)*sizeof(wgttype));

						if ( ret->adjncy == NULL 
							|| ret->adjwgt == NULL )
						{
							printf("Could not allocate");
							printf(" %ld bytes!\n",
								retM_size*sizeof(wgttype));
							abort();
						}
						else
						{
						//	printf("Allocated %d KB!\n",
						//		((retM_size*sizeof(wgttype))/1024));
						}
					}	
/*					if ( i==0 && M_adjncy[l] == 3 )
					{
						printf("Found it, k:%d\n", k);
					}
*/					ret->adjncy[nnzcount-1]=M_adjncy[l];
					ret->adjwgt[nnzcount-1]=(wgttype)wgt*M_adjwgt[l];
					hashtable[M_adjncy[l]]=htcounter++;
				}
				else
				{
					if (ibase+index < 0 || ibase+index>=retM_size) 
					{
						printf("ibase+index:%d\n",ibase+index);
						printf("ibase:%d,i:%d,index:%d\n",
									ibase,i,index);
						abort();
					}	
					wgttype	oldwgt=ret->adjwgt[ibase+index];
					wgttype newwgt;

					newwgt=oldwgt+wgt*M_adjwgt[l];
					ret->adjwgt[ibase+index]=newwgt;
				}
			}
		}
		ret->xadj[i+1]=nnzcount;
		ret->nnz=nnzcount;

		if ( sortLists > 0 )
		{
			ParallelQSort(ret->adjncy, ret->adjwgt,
			ret->xadj[i],ret->xadj[i+1]-1);
		}

		// reset hashtable.
		for(j=ret->xadj[i];j<ret->xadj[i+1];j++)
			hashtable[ret->adjncy[j]]=-1;

/*		if ( i == 250 )
		{
			printf("Yay! Reached 250!\n");
		}
*/
		if ( threshold > 0 )
		{
			int nnzColumn = ret->xadj[i+1] - ret->xadj[i];
			threshPrune( &nnzColumn, (ret->adjncy+ret->xadj[i]),
				(ret->adjwgt+ret->xadj[i]), threshold, 0);
			numPrunedAway += (ret->xadj[i+1]-ret->xadj[i]) -
								nnzColumn;
			nnzcount = ret->xadj[i] + nnzColumn;
			ret->xadj[i+1] = nnzcount;
			ret->nnz = nnzcount;
			if ( nnzColumn == 0 )
				numZeroCols++;
		}

		for ( j = ret->xadj[i]; j<ret->xadj[i+1]; j++)
		{
			if ( ret->adjwgt[j] > max )
				max = ret->adjwgt[j];
			if ( min <= 0 || ret->adjwgt[j] < min )
				min = ret->adjwgt[j];
			if ( threshold > 0 && ret->adjwgt[j] < threshold )
			{
				/*printf("Yikes! ret->adjwgt[%d]:",j);
				printf("%f < threshold: %f, ", ret->adjwgt[j],
				threshold);
				printf("i:%d, adjncy[%d]:%d\n", i, j,
				ret->adjncy[j]); */
			}
			sum += ret->adjwgt[j];
		}

		/*if ( i % 10000 == 0 )
		{
			printf("%d\n",i);
			printf("nnz:%d, numPrunedAway:%ld\n",
			nnzcount, numPrunedAway);
			wgttype avgwgt = sum / nnzcount;
			printf("min:%f, max:%f, avg:%f\n",  min, max, avgwgt);
			printf("numZeroCols:%d\n", numZeroCols);

			fflush(stdout);
		}*/
	 }

	printf("nnzcount of result of expand:%d\n", nnzcount);
	printf("num. pruned away:%ld, threshold:%f\n", 
				numPrunedAway, threshold);
	// dumpMatrix(ret);
	/* for( i=0; i<ret->nvtxs;i++)
		printf("ret_adjwgtsum[i]:%1.5f\n",ret->adjwgtsum[i]);
	 */
	 return ret;

}

void checksums(Matrix* M)
{
	int i,j;
	wgttype sum;
	for(i=0;i<M->nvtxs;i++)
	{
		sum=0;
		for(j=M->xadj[i];j<M->xadj[i+1];j++)
		{
			sum+=M->adjwgt[j];
		}
		if (sum < 0.999 || sum > 1.001 )
		{
			printf("sum not 1. sum:%f, i:%d\n", sum,i);
			abort();
		}
	}
}

/* This method returns the norm of M0-M1. Remember, it assumes
 * that the adjacency lists in both M0 and M1 are sorted in
 * ascending order. */
void changeBetweenMatrices(Matrix* M0, Matrix* M1, wgttype* c)
{
	int i,j,k;
	wgttype sum=0;
//	double sum=0;

	for(i=0;i<M0->nvtxs;i++)
	{
		for(j=M0->xadj[i],k=M1->xadj[i];
			 j<M0->xadj[i+1] && k<M1->xadj[i+1];)
		{
			wgttype a=M0->adjwgt[j],b=M1->adjwgt[k];
			if (M0->adjncy[j]==M1->adjncy[k])
			{
				sum+=(a-b)*(a-b);j++;k++;
			}
			else 
			{
				if (M0->adjncy[j] < M1->adjncy[k])
				{
					sum+=a*a; j++;
				}
				else
				{
					sum+=b*b; k++;
				}
			}
//			printf("sum:%f\n",sum);
		}
		for(;j<M0->xadj[i+1];j++)
			sum+=M0->adjwgt[j]*M0->adjwgt[j];
		for(;k<M1->xadj[i+1];k++)
			sum+=M1->adjwgt[k]*M1->adjwgt[k];

	}

	sum=sqrt(sum);
//	printf("sum:%f\n",sum);
	*c=sum;
//	return sum;
}

/* This function can be used only when the adjacency lists for
 * each vertex are sorted! Beware!
 **/
Matrix* add(Matrix* M0, Matrix* M1)
{
	Matrix* ret;
	int i,j,k; 
	long ret_size=(M0->nnz > M1->nnz)? M0->nnz: M1->nnz ;
	int init_ret_size = ret_size;
	wgttype t;

	ret=allocMatrix(M0->nvtxs, init_ret_size ,0,0,0);
	ret->xadj[0]=0;
	for(i=0;i<M0->nvtxs;i++)
	{
		int l=ret->xadj[i];
		for(j=M0->xadj[i],k=M1->xadj[i];
			 j<M0->xadj[i+1] && k<M1->xadj[i+1];)
		{
			wgttype a=M0->adjwgt[j],b=M1->adjwgt[k];

			if ( l+2 >= ret_size )
			{
				ret_size+=ret_size;
				ret->adjncy=(idxtype*)realloc(ret->adjncy,
								(ret_size)*sizeof(idxtype));
				ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
								(ret_size)*sizeof(wgttype));
				if ( ret->adjncy == NULL || ret->adjwgt == NULL )
				{
					printf("Could not allocate");
					printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
					abort();
				}
			}

			if (M0->adjncy[j]==M1->adjncy[k])
			{
				ret->adjwgt[l]=a+b;
				ret->adjncy[l]=M0->adjncy[j];
				j++;k++;l++;
			}
			else 
			{
				if ( M0->adjncy[j] < M1->adjncy[k] )
				{
					ret->adjwgt[l]=a;
					ret->adjncy[l]=M0->adjncy[j];
					j++;l++;
				}
				else
				{
					ret->adjwgt[l]=b;
					ret->adjncy[l]=M1->adjncy[k];
					k++;l++;
				}
			}
		}
		if ( l+(M0->xadj[i+1]-j)+(M1->xadj[i+1]-k) >= ret_size )	
		{
			ret_size+=ret_size;
			ret->adjncy=(idxtype*)realloc(ret->adjncy,
							(ret_size)*sizeof(idxtype));
			ret->adjwgt=(wgttype*)realloc(ret->adjwgt,
							(ret_size)*sizeof(wgttype));
			if ( ret->adjncy == NULL || ret->adjwgt == NULL )
			{
				printf("Could not allocate");
				printf(" %ld bytes!\n",ret_size*sizeof(idxtype));
				abort();
			}
		}

		for(;j<M0->xadj[i+1];j++)
		{
			ret->adjncy[l]=M0->adjncy[j];
			ret->adjwgt[l++]=M0->adjwgt[j];
		}
		for(;k<M1->xadj[i+1];k++)
		{
			ret->adjncy[l]=M1->adjncy[k];
			ret->adjwgt[l++]=M1->adjwgt[k];
		}

		ret->xadj[i+1]=l;
	}

	ret->nnz=ret->xadj[ret->nvtxs];

	return ret;
}

/*
wgttype flowDifference(int nghbr, Matrix* flows, Matrix* adj,
idxtype* hashtable)
{
	// The hashtable can tell us if our node flows to a
	// particular vertex.
	int i,j;
	wgttype sum=0,k;
	for(i=flows->xadj[nghbr];i<flows->xadj[nghbr+1];i++)
	{
//		if ( flows->adjncy[i] < 0 || flows->adjcn
		j=hashtable[flows->adjncy[i]];
		k=flows->adjwgt[i];
		if ( j > -1 )
		{
			k-=flows->adjwgt[j];
		}
		sum += k*k;
	}
	sum=sqrt(sum);
	return sum;
}
*/

/* The input matrix (flows) to this function needs to have its
 * adjacency lists sorted. This should be in
 * setupCanonicalMatrix. */
wgttype diffTwoColumns(Matrix *flows, int col1, int col2)
{
	int i,j;
	wgttype sum=0;
	for (i=flows->xadj[col1],j=flows->xadj[col2];
			i<flows->xadj[col1+1] && j<flows->xadj[col2+1];)
	{
		wgttype a,b;
		a=flows->adjwgt[i];
		b=flows->adjwgt[j];
		if ( flows->adjncy[i]==flows->adjncy[j] )
		{
			sum += (a-b)*(a-b);
			i++;j++;
		}
		else if ( flows->adjncy[i] < flows->adjncy[j] )
		{
			sum += a*a;
			i++;
		}
		else
		{
			sum += b*b;
			j++;
		}
	}
	for(; i<flows->xadj[col1+1]; i++)
		sum += flows->adjwgt[i]*flows->adjwgt[i];
	for(; j<flows->xadj[col2+1]; j++)
		sum += flows->adjwgt[j]*flows->adjwgt[j];
	
	return sqrt(sum);
}

Matrix* getDprAdjMatrix(Matrix* flows, Matrix* adj, idxtype*
hashtable, wgttype threshold)
{
	int nvtxs=adj->nvtxs,i,j,k,l;
	Matrix* ret=allocMatrix(nvtxs,adj->nnz,0,0,0);
	
/*	for(i=0;i<nvtxs;i++)
		hashtable[i]=-1;
*/
	ret->xadj[0]=0;
	for(i=0;i<nvtxs;i++)
	{
		ret->xadj[i+1]=ret->xadj[i];
		wgttype sum=0;

		// first get the hashtable ready.
/*		for(j=flows->xadj[i];j<flows->xadj[i+1];j++)
		{
			hashtable[flows->adjncy[j]]=j;		
		}
*/
		for(j=adj->xadj[i];j<adj->xadj[i+1];j++)
		{
		/*	if ( flowDifference(adj->adjncy[j],flows,adj,hashtable) 
						< threshold ) */
			if ( diffTwoColumns(flows,i,adj->adjncy[j]) 
						< threshold )
			{
				ret->adjncy[ret->xadj[i+1]]=adj->adjncy[j];
				ret->adjwgt[ret->xadj[i+1]]=adj->adjwgt[j];
				ret->xadj[i+1]++;
				sum+=adj->adjwgt[j];
			}
		}

		// normalize columns
		for(j=ret->xadj[i];j<ret->xadj[i+1];j++)
		{
			ret->adjwgt[j] /= sum;
		}

		// reset the hashtable
/*		for(j=flows->xadj[i];j<flows->xadj[i+1];j++)
			hashtable[flows->adjncy[j]]=-1;
*/		
	}
	ret->nnz=ret->xadj[nvtxs];

	freeMatrix(adj);

	return ret;
}

// Attempt at alleviating many nodes in one cluster problem.
Matrix* transformAdj(Matrix* flows, Matrix* input_adj, 
Options opt)
{
	wgttype* rowSums, *rowSums_2;
	Matrix* adj;
	int i,j;

	rowSums = getRowSums(flows);
	rowSums_2 = (wgttype*)malloc(sizeof(wgttype)*flows->nvtxs);

	for( i=0; i<flows->nvtxs; i++ )
	{
		rowSums_2[i] = 0;
		for ( j=flows->xadj[i]; j<flows->xadj[i+1]; j++ )
		{
			rowSums_2[i] +=	flows->adjwgt[j]
							* rowSums[flows->adjncy[j]];
		}
	}

	free(rowSums);
	rowSums = rowSums_2;

//	printf("penalty_power:%.1f\n", opt.penalty_power);
	adj = allocMatrix(input_adj->nvtxs, input_adj->nnz, 0, 0, 0);
	for ( i=0; i<input_adj->xadj[input_adj->nvtxs]; i++ )
	{
		adj->adjncy[i] = input_adj->adjncy[i];
		wgttype k;
		if ( (k=rowSums[adj->adjncy[i]]) > 1.0 ) 
		{
			if ( opt.penalty_power == 0.5 )
				adj->adjwgt[i] = (1.0/sqrt(k)) * input_adj->adjwgt[i];
			else if ( opt.penalty_power == 1 )
				adj->adjwgt[i] = (1.0/k) * input_adj->adjwgt[i];
			else
			{
				adj->adjwgt[i] = (1.0/pow(k,opt.penalty_power)) *
									input_adj->adjwgt[i];
			}
		}
		else
		{
			adj->adjwgt[i] = input_adj->adjwgt[i];
		}
	}
	idxcopy(adj->nvtxs+1, input_adj->xadj, adj->xadj);
	normalizeColumns(adj, 1, 0); 

	free(rowSums);

	return adj;
}

Matrix* allInOneStep(Matrix* flows, Matrix* input_adj, idxtype*
				hashtable, Options opt,int sortLists, int dprPhase)
{
	int nvtxs=input_adj->nvtxs,i,j,l;	
	long ret_size;
	int init_ret_size;
	int vector_size=100;
	Matrix* ret, *adj;
	idxtype* vector_adj = idxmalloc(vector_size,"Adjacency for a node");
	wgttype* vector_wgt =
			(wgttype*) malloc(sizeof(wgttype)*vector_size);
	float gamma = opt.gamma;
	int exact = opt.exact, k = opt.k;
	
/*	if ( exact <= 0 )
	{
		printf("Don't use allInOneStep without exact pruning!\n");
		abort();
	}
*/	
	for(i=0;i<nvtxs;i++)
	{
		hashtable[i]=-1;
	}

//	printf("In allInOneStep: penalty_power: %.1f\n",
//						opt.penalty_power);
	if ( opt.transformAdj > 0 && opt.penalty_power != 0 && dprPhase < 1 )
		adj=transformAdj(flows, input_adj, opt);
	else
		adj = input_adj;

	//checksums(flows);
//	printf("Successfully checksummed flows\n");
//	checksums(adj);
//	printf("Successfully checksummed adj\n");

	if ( exact > 0 )
		ret_size=nvtxs*k;
	else
		ret_size=flows->nnz;

	init_ret_size = ret_size;
	ret=allocMatrix(nvtxs, init_ret_size,0,0,1);
	ret->xadj[0]=0;
	for(i=0;i<nvtxs;i++)
	{
		int ibase=ret->xadj[i];
		int vector_nnz=0;
		wgttype sum=0,maxwgt=0, initsum,initsum2,thresh;
		
		for(j=adj->xadj[i];j<adj->xadj[i+1];j++)
		{
			wgttype w=adj->adjwgt[j];// weight of neighbour
			int nghbr=adj->adjncy[j];//index of neighbour

			for(l=flows->xadj[nghbr];l<flows->xadj[nghbr+1];l++)
			{
				int new_nghbr=flows->adjncy[l];
				wgttype new_nghbr_wgt=w*flows->adjwgt[l];

				if ( new_nghbr_wgt > 1.01 || w > 1.01 ||
				flows->adjwgt[l] > 1.01 )
				{ 
					printf("oops, wt > 1.");
					printf("new_nghbr_wt:%f",new_nghbr_wgt);
					printf(",w:%f",w);
					printf(",flows->adjwgt:%f\n",flows->adjwgt[l]);
					abort();
				}
				
				if ( new_nghbr < 0 || new_nghbr >= nvtxs )
				{
					printf("Oops! new_nghbr:%d\n",new_nghbr);
					abort();
				}
				int index=hashtable[new_nghbr];
				if ( index == -1 )
				{
					vector_nnz++;
					if ( vector_nnz > vector_size )
					{
						vector_size+=vector_size;
						vector_adj=(idxtype*)realloc(vector_adj,
									sizeof(idxtype)*vector_size);
						vector_wgt=(wgttype*)realloc(vector_wgt,
									sizeof(wgttype)*vector_size);

						if ( vector_wgt == NULL || 
								vector_adj == NULL )
						{
							printf("Could not allocate ");
							printf("%ld KB!\n",
								(vector_size*sizeof(wgttype))/1024);
							abort();
						}
						
					}
					vector_adj[vector_nnz-1]=new_nghbr;
					vector_wgt[vector_nnz-1]=new_nghbr_wgt;
					hashtable[new_nghbr]=vector_nnz-1;
				}
				else
				{
					vector_wgt[index]+=new_nghbr_wgt;
					if ( vector_adj[index] != new_nghbr )
					{
						printf("vector_adj[%d]:%d,",
								index,vector_adj[index]);
						printf(" hashtable says it should be %d\n",
								index);
						printf("vertex id:%d\n",i);
						abort();
					}
					if ( new_nghbr_wgt > 1.01 || vector_wgt[index] >
					1.01 )
					{ 
						printf("oops, wt > 1.");
						printf("new_nghbr_wt:%f",new_nghbr_wgt);
						printf(",vector_wgt:%f\n",vector_wgt[index]);
						abort();
					}	

				}
				
			}
			
		}

		// reset hashtable
		for(j=0;j<vector_nnz;j++)
		{
			if ( vector_adj[j] < 0 || vector_adj[j]>=nvtxs )
			{
				printf("Oops! vector_adj[%d]:%d\n",
							j,vector_adj[j]);
				abort();
			}
			hashtable[vector_adj[j]]=-1;
//			vector_adj[j]=-1;
//			vector_wgt[j]=0;
		}
		
		// check sum before inflating
/*		sum=0;
		for( j=0; j<vector_nnz; j++ )
		{
			sum += vector_wgt[j];
		}
		if ( sum < 0.99 || sum > 1.01 )
		{
			printf("Oops! Sum:%f\n", sum);
			abort();
		}
*/
		// now inflate
		sum=0;
		for(j=0;j<vector_nnz;j++)
		{
			wgttype new_wgt;
			new_wgt=(wgttype)pow((double)vector_wgt[j],(double)gamma);
			vector_wgt[j]=new_wgt;
			sum+=new_wgt;
			
			if ( new_wgt > 1.01 || sum > 1.01 )
			{ 
				printf("oops, wt > 1.");
				printf("new_nghbr_wt:%f",new_wgt);
				printf("sum:%f",sum);
				abort();
			}	

			if ( new_wgt > maxwgt )
				maxwgt=new_wgt;
		}
		 
/*		if ( ! (maxwgt > 0) )
		{
			printf("maxwgt:%f for i:%d\n",maxwgt, i);
			abort();
		}
*/
		thresh=computeThreshold(sum/((wgttype)vector_nnz),maxwgt);
		int old_vector_nnz = vector_nnz;

		initsum=sum;
		initsum2=threshPrune(&vector_nnz,vector_adj,vector_wgt,thresh,sum);
		//I'll directly skip to exact pruning
		if ( exact > 0 )
		{
			sum=exactPrune(&vector_nnz,vector_adj,vector_wgt,k,initsum2);
		}
		else
		{
			sum=initsum2;
		}
		if ( ! (sum > 0)  || sum > 1.01 )
		{
			printf( "Sum:%f for column:%d,initsum:%f",
						sum,i,initsum);
			printf(",initsum2:%f\n",initsum2);
			printf("thresh:%f, maxwgt:%f\n",thresh, maxwgt);
			printf("old_vector_nnz:%d\n", old_vector_nnz);
			abort();
		}

		if ( sortLists > 0 )
		{
			// sort the adjacency lists so that neighbours are
			// listed in increasing order.
			ParallelQSort(vector_adj,vector_wgt,0,vector_nnz-1);
		/*	for(j=0;j<vector_nnz-1;j++)
			{
				if ( vector_adj[j] > vector_adj[j+1] )
				{
					printf("Sort didn't work! ");
					printf("vector_adj[%d]:%d,",j,vector_adj[j]);
					printf("vector_adj[%d]:%d,",j+1,vector_adj[j+1]);
					abort();
				}
			} */
		}

		// now copy vector_adj and normalized vector_wgt to ret
		ret->xadj[i+1]=ret->xadj[i]+vector_nnz;
		ibase=ret->xadj[i];
		ret->attractors[i]=-1;
	
		if ( ret->xadj[i+1] > ret_size )
		{
			ret_size+=ret_size;
			ret->adjncy=(idxtype*)realloc(ret->adjncy,ret_size*sizeof(idxtype));
			ret->adjwgt=(wgttype*)realloc(ret->adjwgt,ret_size*sizeof(wgttype));
			if ( ret->adjncy == NULL 
						|| ret->adjwgt == NULL )
			{
				printf("Could not allocate");
				printf(" %ld bytes!\n",ret_size);
				abort();
			}
		}
		
		/* In MIS-MLRMCL, out-flows of a node can sum to zero,
		 * because if they are not part of the coarser graph,
		 * then their flows are initialized to zero, so for the
		 * first iteration, their out-flows will sum to zero.
		 * Avoid divide by zero error by setting sum to 1. */
		if ( sum == 0 )
			sum=1;

		for(j=ret->xadj[i];j<ret->xadj[i+1];j++)
		{
			ret->adjncy[j]=vector_adj[j-ibase];
			ret->adjwgt[j]=vector_wgt[j-ibase]/sum;
			if ( ret->adjwgt[j] > 0.5 )
				ret->attractors[i]=ret->adjncy[j];
			if ( ret->adjncy[j] < 0 || ret->adjncy[j]>=nvtxs )
			{
				printf("Oops! ret->adjncy[%d]:%d\n",
							j,ret->adjncy[j]);
				abort();
			}
			vector_adj[j-ibase]=-1;
			vector_wgt[j-ibase]=0;
		}

/*		if ( i % 10000 == 0 )
		{
			printf("%d..",i);
			fflush(stdout);
		} */
		
	}
//	printf("\n");

	ret->nnz=ret->xadj[ret->nvtxs];
	free(vector_adj);
	free(vector_wgt);
//	printf("checksumming ret\n");
//	checksums(ret);
//	printf("done checksumming ret\n");

	if ( adj != input_adj )
	{
		freeMatrix(adj);
	}

	return ret;
}

void inflate(Matrix* a, float gamma)
{
	// Will assume that this matrix has adjwgtsum and maxwgt
	// already allocated.
	int i,j,k;
	for(i=0;i<a->nvtxs;i++)
	{
		a->adjwgtsum[i]=0;
		a->maxwgt[i]=0;
		for(j=a->xadj[i];j<a->xadj[i+1];j++)
		{
			wgttype oldwgt, newwgt;
			oldwgt=a->adjwgt[j];
			newwgt=(wgttype)pow((double)oldwgt,(double)gamma);
			a->adjwgt[j]=newwgt;
			a->adjwgtsum[i]+=newwgt;
			if ( newwgt > a->maxwgt[i] )
				a->maxwgt[i]=newwgt;
		}
	}
	
}

void calculateMaxAndSum(Matrix* a)
{
	int i,j;
	if ( a->adjwgtsum == NULL )
		a->adjwgtsum=(wgttype*)malloc(sizeof(wgttype)*a->nvtxs);
	if ( a->maxwgt == NULL )
		a->maxwgt=(wgttype*)malloc(sizeof(wgttype)*a->nvtxs);

	for(i=0;i<a->nvtxs;i++)
	{
		a->adjwgtsum[i]=0;
		a->maxwgt[i]=-1;
		for(j=a->xadj[i];j<a->xadj[i+1];j++)
		{
			a->adjwgtsum[i] += a->adjwgt[j];
			if (a->maxwgt[i]==-1 || a->adjwgt[j] >	a->maxwgt[i])
				a->maxwgt[i] = a->adjwgt[j];
		}
	}
}

void pruneAndNormalize(Matrix* a, int exact, int k)
{
	int i,j,nnz=0,ibase;
	//printf("before prune: avg nnz per column: %d\n",
	//		(a->nnz/a->nvtxs));
//	printf("before prune: nnz:%d, nvtxs:%d\n",a->nnz,a->nvtxs);
	wgttype  *adjwgt, thresh, avg, sum;
	idxtype *xadj, *adjncy,*attractors, attr;
	if ( exact > 0 )
	{
		int toalloc=a->nvtxs*k;
		toalloc=(a->nnz < toalloc) ? a->nnz : toalloc;
		adjncy=idxmalloc(toalloc,"Adjacency data");
		adjwgt=(wgttype*)malloc(sizeof(wgttype)*toalloc);
	}
	else
	{
		adjncy=idxmalloc(a->nnz,"Adjacency data");
		adjwgt=(wgttype*)malloc(sizeof(wgttype)*(a->nnz));
	}
	xadj=idxmalloc(a->nvtxs+1,"Xadj data");
	attractors=idxmalloc(a->nvtxs,"Attractors");

	if ( a->adjwgtsum == NULL || a->maxwgt == NULL )
		calculateMaxAndSum(a);

	xadj[0]=0;
//	for( i=0; i<a->nnz;i++)
//		printf("M_adjwgt[%d]:%1.5f\n",i,a->adjwgt[i]);

/*	for(i=0;i<a->nvtxs;i++)
	{
		fprintf(stderr,"a->maxwgt[%d]:%1.4f,",i,a->maxwgt[i]);
		fprintf(stderr,"a->wgtsum[%d]:%1.4f\n",i,a->adjwgtsum[i]);
	}
*/
	for (i=0; i<a->nvtxs; i++)
	{
		int degree;
		wgttype maxwgt;
		ibase=xadj[i];
		degree=a->xadj[i+1]-a->xadj[i];
		if ( degree > 0 )
			avg=a->adjwgtsum[i]/degree;
		else
			continue;

		thresh=computeThreshold(avg, a->maxwgt[i]);
		sum=0;
		for(j=a->xadj[i];j<a->xadj[i+1];j++)
		{
			if ( a->adjwgt[j]>=thresh )
			{
				adjncy[nnz]=a->adjncy[j];
				adjwgt[nnz]=a->adjwgt[j];
				nnz++;
				sum+=a->adjwgt[j];
			}
		}	
		xadj[i+1]=nnz;

		if ( exact > 0 )
		{
			int nnghbr=xadj[i+1]-xadj[i];
			wgttype initsum=sum;
			sum = exactPrune(&nnghbr, (adjncy+xadj[i]),
					(adjwgt+xadj[i]), k, sum);
			xadj[i+1]=nnghbr+xadj[i];
			nnz=xadj[i+1];
			if ( ! (sum > 0) )
			{
				fprintf(stderr, "Sum:%f for column:%d,",sum,i);
				fprintf(stderr, "k:%d, nnghbr:%d,initsum:%f\n",k,
						nnghbr, initsum);
	/*			for(j=a->xadj[i];j<a->xadj[i+1];j++)
					fprintf(stderr,"a->adjwgt[%d]:%1.5f\n",a->adjncy[j],a->adjwgt[j]);*/
				abort();
			}
		}

		if ( ! (sum > 0) )
		{
			fprintf(stderr, "Sum:%f for column:%d,",sum,i);
			fprintf(stderr, "thresh:%f\n", thresh);
/*			for(j=a->xadj[i];j<a->xadj[i+1];j++)
				fprintf(stderr,"a->adjwgt[%d]:%1.5f\n",a->adjncy[j],a->adjwgt[j]);*/
			abort();
		}

		attr=-1;
		maxwgt=0;
		for(j=ibase;j<nnz;j++)
		{
			adjwgt[j] /= sum;
			if ( adjwgt[j] > 0.5 )
			//if ( adjwgt[j] > maxwgt )
			{
				attr = adjncy[j];;
			//	maxwgt = adjwgt[j];
			}
		}
		attractors[i]=attr;
		
	}

	free(a->adjncy);
	free(a->adjwgt);
	free(a->xadj);
	free(a->adjwgtsum); /* going to be 1 anyway, don't need this
	*/
	free(a->maxwgt); /* don't need this during expand */
	if ( a->attractors != NULL )
		free(a->attractors);
		
	a->adjwgt=adjwgt;
	a->adjncy=adjncy;
	/*a->adjwgt=realloc(adjwgt,nnz);
	a->adjncy=realloc(adjncy,nnz); */
	a->attractors=attractors;
	a->nnz=nnz;
	a->xadj=xadj;
	a->adjwgtsum=NULL;
	a->maxwgt=NULL;

	//printf("after prune: avg nnz per column: %d\n",
	//		(a->nnz/a->nvtxs));
//	printf("after prune, nnz:%d\n",nnz);

}


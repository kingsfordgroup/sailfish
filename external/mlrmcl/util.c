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


/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * util.c
 *
 * This function contains various utility routines
 *
 * Started 9/28/95
 * George
 *
 * Significant code has been added by Venu Satuluri
 *
 * $Id: util.c,v 1.17 2010-11-30 23:21:46 venu Exp $
 */

#include <metis.h>

void initOptions(Options* opt)
{
	opt->gamma = 2.0;
	opt->coarsenTo = 1000;
	opt->iter_per_level = 4;
	opt->num_last_iter = 25;
	opt->ncutify = 1;
	opt->hubRemoval = 0;
	opt->hubPct = 0.01;
	opt->mis_coarsenType = 0;
	opt->exact = 0;
	opt->k = 100;
	opt->dpr_threshold = 0.4;
	opt->transformAdj = 1;
//	opt->matchType = MATCH_SHEMN;
	opt->matchType = MATCH_UNSPECIFIED;
	opt->penalty_power = 0.5;
}

void mapIndices(idxtype *indices, idxtype* map, int n, int npart)
{
	int i=n-1;
	int extraIndexCounter=0;
	for ( ; i >= 0; i-- )
	{
		if ( map[i] > -1 )
			indices[i] = indices[map[i]];
		else
			indices[i] = npart + (extraIndexCounter++);
	}

}

void assignClustersToHubs(idxtype* indices, idxtype* map, int n,
int npart, GraphType* graph)
{
	idxtype* clusterCount = idxmalloc(npart,
						"assignClustersToHubs:clusterCount");
	idxtype* newIndices = idxmalloc(graph->nvtxs,
						"assignClustersToHubs:newIndices");

	for ( int i=0; i<n; i++ )
	{
		if ( map[i] > -1 )
		{
			newIndices[i] = indices[map[i]];
			continue;
		}
		for ( int j=0; j<npart; j++ )
			clusterCount[j] = 0;
		for ( int k=graph->xadj[i]; k<graph->xadj[i+1]; k++ )
		{
			if ( map[graph->adjncy[k]] > -1 )
				clusterCount[indices[map[graph->adjncy[k]]]]++;
		}
		int bestCluster = 0, bestClusterCount = 0;
		for ( int j=0; j<npart; j++ )
		{
			if ( clusterCount[j] > bestClusterCount )
			{
				bestCluster = j;
				bestClusterCount = clusterCount[j];
			}
		}
		newIndices[i] = bestCluster;
	}

	for ( int i=0; i<graph->nvtxs; i++ )
		indices[i] = newIndices[i];

	free(newIndices);
}


idxtype* lookForSingletons(GraphType* graph, int* noOfSingletons)
{
	idxtype* newIds = idxmalloc(graph->nvtxs, 
							"lookForSingletons:newIds");
	int i, newIdCounter=0;

	for ( i=0; i<graph->nvtxs; i++ )
	{
		newIds[i]=newIdCounter++;
		if ( graph->xadj[i+1]==graph->xadj[i] 
		 || ( graph->xadj[i+1]==graph->xadj[i]+1 
		   && graph->adjncy[graph->xadj[i]] == i ) )
		{
			newIds[i] = -1;
			newIdCounter--;
		}
	}

	*noOfSingletons = graph->nvtxs - newIdCounter;

	if ( *noOfSingletons > 0 )
		return newIds;
	else
	{
		free(newIds);
		return NULL;
	}

}

void dfTraversalMatrix(Matrix *graph, idxtype root, idxtype*
visited, int* nVisited, wgttype minWgt)
{
	int i;
	if ( visited[root] < 1 )
	{
		visited[root] = 1;
		(*nVisited)++;
	}
	for ( i=graph->xadj[root]; i<graph->xadj[root+1]; i++ )
	{
		 if ( graph->adjwgt[i] >= minWgt && 
		 		visited[graph->adjncy[i]] < 1 )
		 {
		 	visited[graph->adjncy[i]]=1;
			(*nVisited)++;
			dfTraversalMatrix(graph, graph->adjncy[i], visited,
			nVisited, minWgt);
		 }
	}
}


void dfTraversal(GraphType *graph, idxtype root, idxtype* visited, int* nVisited)
{
	int i;
	if ( visited[root] < 1 )
	{
		visited[root] = 1;
		(*nVisited)++;
	}
	for ( i=graph->xadj[root]; i<graph->xadj[root+1]; i++ )
	{
		 if ( visited[graph->adjncy[i]] < 1 )
		 {
		 	visited[graph->adjncy[i]]=1;
			(*nVisited)++;
			dfTraversal(graph, graph->adjncy[i], visited,
			nVisited);
		 }
	}
}

/*************************************************************************
* This function prints an error message and exits
**************************************************************************/
void errexit(char *f_str,...)
{
  va_list argp;
  char out1[256], out2[256];

  va_start(argp, f_str);
  vsprintf(out1, f_str, argp);
  va_end(argp);

//  sprintf(out2, "Error! %s", out1);

  printf("Error! %s", out1);
  fflush(stdout);

  abort();
}



#ifndef DMALLOC
/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
int *imalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (int *)GKmalloc(sizeof(int)*n, msg);
}


/*************************************************************************
* The following function allocates an array of integers
**************************************************************************/
idxtype *idxmalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (idxtype *)GKmalloc(sizeof(idxtype)*n, msg);
}


/*************************************************************************
* The following function allocates an array of float 
**************************************************************************/
float *fmalloc(int n, const char *msg)
{
  if (n == 0)
    return NULL;

  return (float *)GKmalloc(sizeof(float)*n, msg);
}


/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
int *ismalloc(int n, int ival, const char *msg)
{
  if (n == 0)
    return NULL;

  return iset(n, ival, (int *)GKmalloc(sizeof(int)*n, msg));
}



/*************************************************************************
* The follwoing function allocates an array of integers
**************************************************************************/
idxtype *idxsmalloc(int n, idxtype ival, const char *msg)
{
  if (n == 0)
    return NULL;

  return idxset(n, ival, (idxtype *)GKmalloc(sizeof(idxtype)*n, msg));
}

/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(long nbytes, const char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
  {
    printf("In GKmalloc for long\n");
    errexit("***Memory allocation failed for %s. Requested size: %ld bytes", msg, nbytes);
  }

  return ptr;
}

/*************************************************************************
* This function is my wrapper around malloc
**************************************************************************/
void *GKmalloc(int nbytes, const char *msg)
{
  void *ptr;

  if (nbytes == 0)
    return NULL;

  ptr = (void *)malloc(nbytes);
  if (ptr == NULL) 
  {
    printf("In GKmalloc for int\n");
    errexit("***Memory allocation failed for %s. Requested size: %d bytes", msg, nbytes);
  }
  
  return ptr;
}
#endif

/*************************************************************************
* This function is my wrapper around free, allows multiple pointers    
**************************************************************************/
void GKfree(void **ptr1,...)
{
  va_list plist;
  void **ptr;

  if (*ptr1 != NULL)
    free(*ptr1);
  *ptr1 = NULL;

  va_start(plist, ptr1);

  /* while ((int)(ptr = va_arg(plist, void **)) != -1) { */
  while ((ptr = va_arg(plist, void **)) != LTERM) {
    if (*ptr != NULL)
      free(*ptr);
    *ptr = NULL;
  }

  va_end(plist);
}            


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
int *iset(int n, int val, int *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
idxtype *idxset(int n, idxtype val, idxtype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}

wgttype* wgtset(int n, wgttype val, wgttype *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}


/*************************************************************************
* These functions set the values of a vector
**************************************************************************/
float *sset(int n, float val, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] = val;

  return x;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int iamax(int n, int *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}


/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax(int n, idxtype *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int idxamax_strd(int n, idxtype *x, int incx)
{
  int i, max=0;

  n *= incx;
  for (i=incx; i<n; i+=incx)
    max = (x[i] > x[max] ? i : max);

  return max/incx;
}



/*************************************************************************
* These functions return the index of the maximum element in a vector
**************************************************************************/
int samax(int n, float *x)
{
  int i, max=0;

  for (i=1; i<n; i++)
    max = (x[i] > x[max] ? i : max);

  return max;
}

/*************************************************************************
* These functions return the index of the almost maximum element in a vector
**************************************************************************/
int samax2(int n, float *x)
{
  int i, max1, max2;

  if (x[0] > x[1]) {
    max1 = 0;
    max2 = 1;
  }
  else {
    max1 = 1;
    max2 = 0;
  }

  for (i=2; i<n; i++) {
    if (x[i] > x[max1]) {
      max2 = max1;
      max1 = i;
    }
    else if (x[i] > x[max2])
      max2 = i;
  }

  return max2;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int idxamin(int n, idxtype *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* These functions return the index of the minimum element in a vector
**************************************************************************/
int samin(int n, float *x)
{
  int i, min=0;

  for (i=1; i<n; i++)
    min = (x[i] < x[min] ? i : min);

  return min;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum(int n, idxtype *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int idxsum_strd(int n, idxtype *x, int incx)
{
  int i, sum = 0;

  for (i=0; i<n; i++, x+=incx) {
    sum += *x;
  }

  return sum;
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void idxadd(int n, idxtype *x, idxtype *y)
{
  for (n--; n>=0; n--)
    y[n] += x[n];
}


/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int charsum(int n, const char *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
int isum(int n, int *x)
{
  int i, sum = 0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum(int n, float *x)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++)
    sum += x[i];

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
float ssum_strd(int n, float *x, int incx)
{
  int i;
  float sum = 0.0;

  for (i=0; i<n; i++, x+=incx)
    sum += *x;

  return sum;
}

/*************************************************************************
* This function sums the entries in an array
**************************************************************************/
void sscale(int n, float alpha, float *x)
{
  int i;

  for (i=0; i<n; i++)
    x[i] *= alpha;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float snorm2(int n, float *v)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += v[i] * v[i];

  return sqrt(partial);
}



/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
float sdot(int n, float *x, float *y)
{
  int i;
  float partial = 0;
 
  for (i = 0; i<n; i++)
    partial += x[i] * y[i];

  return partial;
}


/*************************************************************************
* This function computes a 2-norm
**************************************************************************/
void saxpy(int n, float alpha, float *x, int incx, float *y, int incy)
{
  int i;
 
  for (i=0; i<n; i++, x+=incx, y+=incy) 
    *y += alpha*(*x);
}



/** Venu: my addition here */
/* This function returns the vertices of the graph sorted in
 * ascending order according to the number of neighbours they
 * have. The idea is that during matching, if the vertices with
 * low degrees are matched first, then there is lesser risk of
 * some vertices being "orphaned". */
void permuteDegreeOrder(int n, idxtype *p, idxtype *xadj)
{
	idxtype* degrees=idxmalloc(n,"permuteDegreeOrder:degrees");
	wgttype* temp_p=(wgttype*)malloc(sizeof(wgttype)*n);
	int i;

	for( i=0; i<n; i++ )
	{
		degrees[i]=xadj[i+1]-xadj[i];
		temp_p[i]=i;
	}

	ParallelQSort(degrees,temp_p,0,n-1);

	for ( i=0; i<n; i++ )
		p[i]=(idxtype)temp_p[i];

	free(degrees);
	free(temp_p);
}


/*************************************************************************
* This file randomly permutes the contents of an array.
* flag == 0, don't initialize perm
* flag == 1, set p[i] = i 
**************************************************************************/
void RandomPermute(int n, idxtype *p, int flag)
{
  int i, u, v;
  idxtype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for ( i=0; i<n; i++ )
  {
  	u = RandomInRangeFast(n);
	SWAP(p[i], p[u], tmp);
  }

/*  for (i=0; i<n; i+=16) {
    u = RandomInRangeFast(n-4);
    v = RandomInRangeFast(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  } */
}

void RandomPermuteWgttype(int n, wgttype *p, int flag)
{
  int i, u, v;
  wgttype tmp;

  if (flag == 1) {
    for (i=0; i<n; i++)
      p[i] = i;
  }

  if (n <= 4)
    return;

  for (i=0; i<n; i+=16) {
    u = RandomInRangeFast(n-4);
    v = RandomInRangeFast(n-4);
    SWAP(p[v], p[u], tmp);
    SWAP(p[v+1], p[u+1], tmp);
    SWAP(p[v+2], p[u+2], tmp);
    SWAP(p[v+3], p[u+3], tmp);
  }
}

/* The specification for this function is a bit complex. The
 * input array a is sorted; if 'key' is present in a, then this
 * function is supposed to work like normal binary search and
 * return the position of 'key' in a. Otherwise, it is supposed to
 * return the position at which key should be inserted if the
 * array a is to remain sorted. */
int bsearch_insertPos(idxtype *a, idxtype key, int start, int end)
{
	if ( start > end )
		return start;

	if ( key == a[start] )
		return start;
	if ( key == a[end] )
		return end;
	if ( end - start < 2 )
	{
		/* These 3 cases should cover all possiblities. */
		if ( key < a[start] )
			return start;
		if ( key < a[end] )
			return end;
		if ( key > a[end] )
			return end+1;
	}

	int mid = (start+end)/2;
	if ( mid+1 > end || mid < start )
	{
		printf("start:%d, mid:%d, end:%d\n", start, mid, end);
		abort();
	}
	if ( key == a[mid] )
		return mid;
	if ( key < a[mid] )
		return bsearch_insertPos(a, key, start, mid);
	else
		return bsearch_insertPos(a, key, mid+1, end);

}


int RandomPartitionInts(idxtype* a, int start, int end)
{
	int n=end-start+1,i,j;
	int tmp,x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
		SWAP(a[i],a[end],tmp);

	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmp);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmp);

	return i+1;
}

int RandomPartition(wgttype* a, int start, int end)
{
	int n=end-start+1,i,j;
	wgttype tmp,x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
		SWAP(a[i],a[end],tmp);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmp);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmp);

	return i+1;
}

int ParallelRandomPartitionInts(idxtype* a, idxtype* b, int start,
									int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

void ParallelQSortInts(idxtype *a, idxtype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionInts(a,b,start,end);
		ParallelQSortInts(a,b,start,q-1);
		ParallelQSortInts(a,b,q+1,end);
	}
}

int ParallelRandomPartitionIntsUsingScores(idxtype* a, idxtype*
				b, idxtype *scores,	int start, int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

	//x=a[end];
	x = scores[a[end]];
	i=start-1;
	for(j=start;j<end;j++)
	{
		//if (a[j] <= x)
		if ( scores[a[j]] <= x )
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

// we will use scores as the comparer to sort a and b together. 
void ParallelQSortIntsUsingScores(idxtype *a, idxtype *b, idxtype
*scores, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionIntsUsingScores(a,b,scores,start,end);
		ParallelQSortIntsUsingScores(a,b,scores,start,q-1);
		ParallelQSortIntsUsingScores(a,b,scores,q+1,end);
	}
}

int RandomPartitionIntsUsingInts(idxtype* a, const idxtype *scores,
				int start, int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	idxtype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
	}

	//x=a[end];
	x = scores[a[end]];
	i=start-1;
	for(j=start;j<end;j++)
	{
		//if (a[j] <= x)
		if ( scores[a[j]] <= x )
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);

	return i+1;
}

// we will use scores as the comparer to sort a. 
void QSortIntsUsingInts(idxtype *a, idxtype *scores, int start,
int end)
{
	int q;
	if ( start < end )
	{
		q=RandomPartitionIntsUsingInts(a,scores,start,end);
		QSortIntsUsingInts(a,scores,start,q-1);
		QSortIntsUsingInts(a,scores,q+1,end);
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartitionFloatsInts(wgttype* a, idxtype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	int tmpidx;
	wgttype tmpwgt, x;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	

	if ( i != end )
	{
		SWAP(a[i],a[end],tmpwgt);
		SWAP(b[i],b[end],tmpidx);
	}

	x=a[end];
	i=start-1;
	for(j=start;j<end;j++)
	{
		if (a[j] <= x)
		{
			i++;
			if ( i == j )
				continue;
			SWAP(a[i],a[j],tmpwgt);
			SWAP(b[i],b[j],tmpidx);
		}
	}
	SWAP(a[i+1],a[end],tmpwgt);
	SWAP(b[i+1],b[end],tmpidx);

	return i+1;
}

/* sorts both a and b, according to the order imposed by a.*/
void ParallelQSortFloatsInts(wgttype *a, idxtype *b, int start,
int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionFloatsInts(a,b,start,end);
		ParallelQSortFloatsInts(a,b,start,q-1);
		ParallelQSortFloatsInts(a,b,q+1,end);
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartitionLongs(long* a, wgttype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	long tmpidx,x;
	wgttype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

/* Sorts the arrays a and b, according to the order imposed by a.
   [start,end] is the inclusive range of the two arrays i.e.
   the length of the two arrays is end-start+1, not end-start. */
void ParallelQSortLongs(long *a, wgttype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartitionLongs(a,b,start,end);
		ParallelQSortLongs(a,b,start,q-1);
		ParallelQSortLongs(a,b,q+1,end);
	}
}

/* This routine partitions both the arrays a and b, according to
 * the order imposed by a. */
int ParallelRandomPartition(idxtype* a, wgttype* b, int start,
int end)
{
	int n=end-start+1,i,j;
	idxtype tmpidx,x;
	wgttype tmpwgt;
	if ( n<2)
		return start;
	i=start+RandomInRangeFast(n);	
//	printf("chose i:%d, end:%d\n",i,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	if ( i != end )
	{
		SWAP(a[i],a[end],tmpidx);
		SWAP(b[i],b[end],tmpwgt);
	}

/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
*/
	x=a[end];
//	printf("x:%f\n",x);
	i=start-1;
	for(j=start;j<end;j++)
	{
//		printf("Got here, %f\n");
		if (a[j] <= x)
		{
			i++;
			SWAP(a[i],a[j],tmpidx);
			SWAP(b[i],b[j],tmpwgt);
//			printf("i:%d, a[i]:%f\n",i,a[i]);
		}
	}
	SWAP(a[i+1],a[end],tmpidx);
	SWAP(b[i+1],b[end],tmpwgt);

	return i+1;
}

/* Sorts the arrays a and b, according to the order imposed by a.
   [start,end] is the inclusive range of the two arrays i.e.
   the length of the two arrays is end-start+1, not end-start. */
void ParallelQSort(idxtype *a, wgttype *b, int start, int end)
{
	int q;
	if ( start < end )
	{
		q=ParallelRandomPartition(a,b,start,end);
		ParallelQSort(a,b,start,q-1);
		ParallelQSort(a,b,q+1,end);
	}
}

idxtype RandomSelectInts(idxtype *a, int start, int end, int i)
{
	int q,k,j;
	if (start==end)
		return a[start];
	
	q = RandomPartitionInts(a,start,end);

//	printf("After rand partition, q:%d, end:%d\n",q,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
	abort();
*/
	k=q-start+1;
	if (k==i)
		return a[q];
	if (i < k )
		return RandomSelectInts(a,start,q-1,i);
	else
		return RandomSelectInts(a,q+1,end,i-k);
}

wgttype RandomSelect(wgttype *a, int start, int end, int i)
{
	int q,k,j;
	if (start==end)
		return a[start];
	
	q=RandomPartition(a,start,end);

//	printf("After rand partition, q:%d, end:%d\n",q,end);
/*	for(j=start;j<=end;j++)
		printf("%f,",a[j]);
	printf("\n");
	abort();
*/
	k=q-start+1;
	if (k==i)
		return a[q];
	if (i < k )
		return RandomSelect(a,start,q-1,i);
	else
		return RandomSelect(a,q+1,end,i-k);
}


/*************************************************************************
* This function returns true if the a is a power of 2
**************************************************************************/
int ispow2(int a)
{
  for (; a%2 != 1; a = a>>1);
  return (a > 1 ? 0 : 1);
}


/*************************************************************************
* This function initializes the random number generator
**************************************************************************/
void InitRandom(int seed)
{
  if (seed == -1) {
#ifndef __VC__
    srand48(7654321L);  
#endif
    srand(4321);  
  }
  else {
#ifndef __VC__
    srand48(seed);  
#endif
    srand(seed);  
  }
}

/*************************************************************************
* This function returns the log2(x)
**************************************************************************/
int log2(int a)
{
  int i;

  for (i=1; a > 1; i++, a = a>>1);
  return i-1;
}

/*************************************************************************
* This function clears the timers
**************************************************************************/
void InitTimers(CtrlType *ctrl)
{
  cleartimer(ctrl->TotalTmr);
  cleartimer(ctrl->InitPartTmr);
  cleartimer(ctrl->MatchTmr);
  cleartimer(ctrl->ContractTmr);
  cleartimer(ctrl->CoarsenTmr);
  cleartimer(ctrl->UncoarsenTmr);
  cleartimer(ctrl->RefTmr);
  cleartimer(ctrl->ProjectTmr);
  cleartimer(ctrl->SplitTmr);
  cleartimer(ctrl->SepTmr);
  cleartimer(ctrl->AuxTmr1);
  cleartimer(ctrl->AuxTmr2);
  cleartimer(ctrl->AuxTmr3);
  cleartimer(ctrl->AuxTmr4);
  cleartimer(ctrl->AuxTmr5);
  cleartimer(ctrl->AuxTmr6);
}



/*************************************************************************
* This function prints the various timers
**************************************************************************/
void PrintTimers(CtrlType *ctrl)
{
  printf("\nTiming Information -------------------------------------------------");
  printf("\n Multilevel: \t\t %7.3f", gettimer(ctrl->TotalTmr));
  printf("\n     Coarsening: \t\t %7.3f", gettimer(ctrl->CoarsenTmr));
  printf("\n            Matching: \t\t\t %7.3f", gettimer(ctrl->MatchTmr));
  printf("\n            Contract: \t\t\t %7.3f", gettimer(ctrl->ContractTmr));
  printf("\n     Initial Partition: \t %7.3f", gettimer(ctrl->InitPartTmr));
  printf("\n   Construct Separator: \t %7.3f", gettimer(ctrl->SepTmr));
  printf("\n     Uncoarsening: \t\t %7.3f", gettimer(ctrl->UncoarsenTmr));
  printf("\n          Refinement: \t\t\t %7.3f", gettimer(ctrl->RefTmr));
  printf("\n          Projection: \t\t\t %7.3f", gettimer(ctrl->ProjectTmr));
  printf("\n     Splitting: \t\t %7.3f", gettimer(ctrl->SplitTmr));
  printf("\n          AUX1: \t\t %7.3f", gettimer(ctrl->AuxTmr1));
  printf("\n          AUX2: \t\t %7.3f", gettimer(ctrl->AuxTmr2));
  printf("\n          AUX3: \t\t %7.3f", gettimer(ctrl->AuxTmr3));
  printf("\n********************************************************************\n");
}


/*************************************************************************
* This function returns the seconds
**************************************************************************/
double seconds(void)
{
  return((double) clock()/CLOCKS_PER_SEC);
}

/*************************************************************************
* This function uses simple counting sort to return a permutation array
* corresponding to the sorted order. The keys are assumed to start from
* 0 and they are positive.  This sorting is used during matching.
**************************************************************************/
void BucketSortKeysInc(int n, int max, idxtype *keys, idxtype *tperm, idxtype *perm)
{
  int i, ii;
  idxtype *counts;

  counts = idxsmalloc(max+2, 0, "BucketSortKeysInc: counts");

  for (i=0; i<n; i++)
    counts[keys[i]]++;
  MAKECSR(i, max+1, counts);

  for (ii=0; ii<n; ii++) {
    i = tperm[ii];
    perm[counts[keys[i]]++] = i;
  }

  free(counts);
}

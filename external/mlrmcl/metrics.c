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
 *
 * $Id: debug.c,v 1.18 2010/02/04 04:46:27 venu Exp $
 *
 */

#include <metis.h>

int ComputeNumEdgesCut(GraphType* graph, idxtype* where)
{
  	int i, j, cut;
	idxtype* nodeCuts = idxmalloc(graph->nvtxs,
								"ComputeNumEdgesCut:nodeCuts");

	for ( i=0; i<graph->nvtxs; i++ )
		nodeCuts[i] = 0;

    for (cut=0, i=0; i<graph->nvtxs; i++) 
	{
    	for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
		{
        	if ( where[i] != where[graph->adjncy[j]] )
			{
          		cut++;
				if ( where[i] > 0 )
					nodeCuts[i]++;
	/*			if ( (where[i] != 1 && where[graph->adjncy[j]]!= 1)
				 || (where[i] != 0 && where[graph->adjncy[j]] !=
				 0) )
				{
					printf("%d %d\n",
							where[i], where[graph->adjncy[j]]);
					abort(); 
				}*/
			}
		}
    }

	int testCut=0;
	for ( i=0; i<graph->nvtxs; i++ )
		testCut += nodeCuts[i];

//	printf("1. testCut:%d, cut:%d\n", testCut, cut/2);

	free(nodeCuts);
 
 	return cut;
}


/*************************************************************************
* This function computes the cut given the graph and a where vector
**************************************************************************/
int ComputeCut(GraphType *graph, idxtype *where)
{
  int i, j, cut;

  if (graph->adjwgt == NULL) {
//  printf("adjwgt null\n");
  	cut = ComputeNumEdgesCut(graph, where);
   }
  else {
    for (cut=0, i=0; i<graph->nvtxs; i++) {
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        if (where[i] != where[graph->adjncy[j]])
          cut += graph->adjwgt[j];
    }
  }

  return cut/2;
}

void ComputeAdjWgtSums(GraphType* graph)
{
	if ( graph->adjwgtsum == NULL )
	{
		graph->adjwgtsum = idxmalloc(graph->nvtxs,
			"ComputeAdjWgtSums:graph->adjwgtsum");
	}

	idxtype sum;
	int i, j;
    for (i=0; i<graph->nvtxs; i++) 
	{
      sum = 0;
      for (j=graph->xadj[i]; j<graph->xadj[i+1]; j++)
        sum += (graph->adjwgt == NULL) ? 1 : graph->adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }
}

/* compute volume of subset of graph for which bisection[i]>0 */
idxtype ComputeVolSubset(GraphType* graph, idxtype* bisection)
{
	int i,j;
	idxtype vol=0; 
	for ( i=0; i<graph->nvtxs; i++ )
	{
		if ( bisection[i] > 0 )
		{
	/*		if ( graph->adjwgtsum != NULL )
				vol += graph->adjwgtsum[i];
			else
			{
	*/			if ( graph->adjwgt != NULL )
				{
					int k=0;
					for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
					{
						if ( graph->adjncy[j] != i )
							k += graph->adjwgt[j];		
					}
					vol += k;
				}
				else
				{
					int k=0;
					for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
					{
						if ( graph->adjncy[j] != i )
							k++;
					}
					vol += k;
				}
	//		}
		}
	}

	return vol;
}

/*
void computeDirCutAndDirVol(Matrix *graph, wgttype* pagerank, 
	idxtype* members, int size, idxtype* ht, wgttype* cut, wgttype* vol)
{
	int i,j;
	*cut = 0;
	*vol = 0;

	for ( i = 0; i<size; i++ )
	{
		for ( j = graph->xadj[members[i]];
				j<graph->xadj[members[i]+1]; j++ )
		{
			if ( graph->adjncy[j] != members[i] )
			{
				(*vol) += pagerank[members[i]]*graph->adjwgt[j];
				if ( ht[graph->adjncy[j]] <= 0 )
					(*cut) += pagerank[members[i]]*graph->adjwgt[j];
			}
		}
	}
	
}
*/

void computeCutAndVol(GraphType* graph, idxtype* members, int
size, idxtype* ht, int* cut, int *vol)
{
	int i,j;
	*cut = 0;
	*vol = 0;
	if ( graph->adjwgt != NULL )
	{
		for ( i = 0; i<size; i++ )
		{
			for ( j = graph->xadj[members[i]];
					j<graph->xadj[members[i]+1]; j++ )
			{
				if ( graph->adjncy[j] != members[i] )
				{
					(*vol) += graph->adjwgt[j];
					if ( ht[graph->adjncy[j]] <= 0 )
						(*cut) += graph->adjwgt[j];
				}
			}
		}
	}
	else
	{
		for ( i = 0; i<size; i++ )
		{
			for ( j = graph->xadj[members[i]];
					j<graph->xadj[members[i]+1]; j++ )
			{
				if ( graph->adjncy[j] != members[i] )
				{
					(*vol)++;
					if ( ht[graph->adjncy[j]] <= 0 )
						(*cut)++;
				}
			}
		}
	}
}

/*
float Overlap_ComputeNcutVector(GraphType *graph, idxtype
	*members, idxtype *groups_xadj, idxtype ngroups, float *ncutVector)
{
	int i, j;
	idxtype *ht = idxmalloc(graph->nvtxs,
							"Overlap_ComputeNcutVector:ht");
	long graphVol = 0;
	wgttype dirGraphVol = 0;
	Matrix *M;

	for ( i=0; i<graph->nvtxs; i++ )
		ht[i]=0;

	if ( graph->isDirected > 0 )
	{
		if ( graph->pagerank == NULL )
		{
			PageRankOptions opt;
			initPageRankOptions(&opt);
			graph->pagerank = pagerank(graph, opt);
		}

		M = setupCanonicalMatrix( graph->nvtxs, graph->nedges,
		graph->xadj, graph->adjncy, graph->adjwgt, 0); // last arg (0)
		// specifies no ncutify stuff.	

		// if pagerank is computed correctly, dirGraphVol is
		//  simply sum of pageranks i.e. 1. 
		 dirGraphVol = 1.0;
	}
	else
	{
		if ( graph->adjwgt != NULL )
		{
			for ( i=0; i<graph->xadj[graph->nvtxs]; i++ )
				graphVol += graph->adjwgt[i];
		}
		else
			graphVol = graph->xadj[graph->nvtxs];
	}

	float avgCut = 0;
	for ( i=0; i<ngroups; i++ )
	{
		for ( j=groups_xadj[i]; j < groups_xadj[i+1]; j++	)
			ht[members[j]] = 1;

//		int cut = ComputeCut( graph, ht);
//		int vol = ComputeVolSubset( graph, ht);

		int cut, vol;
		wgttype dirCut, dirVol;

		if ( graph->isDirected > 0 )
		{
			computeDirCutAndDirVol(M, graph->pagerank, 
			members+groups_xadj[i], 
			groups_xadj[i+1]-groups_xadj[i], ht, &dirCut, &dirVol); 

			dirVol = (dirVol > (dirGraphVol - dirVol)) ? (dirGraphVol-dirVol) : dirVol;
			ncutVector[i] = (dirVol == 0) ? 0: (dirCut) / (dirVol);
		}
		else
		{
			computeCutAndVol( graph, members+groups_xadj[i],
			groups_xadj[i+1]-groups_xadj[i], ht, &cut, &vol);

			vol = ( vol > (graphVol-vol) )? (graphVol-vol) : vol ;

			ncutVector[i] = (vol==0) ? 0 : ((float) cut) / ((float) vol);
		}

		avgCut += (ncutVector[i] / (float)ngroups);

		for ( j=groups_xadj[i]; j < groups_xadj[i+1]; j++	)
			ht[members[j]] = 0;

	}

//	avgCut /= ((float)ngroups);

	GKfree( (void**)&ht, LTERM);
	return avgCut;
}
*/

/*************************************************************************
* This function computes the normalized cut given the graph and a where vector
**************************************************************************/
float ComputeNCutVector(GraphType *graph, idxtype *where, int
npart,float* ncutVector)
{
  int i, j, cm, nvtxs;
  idxtype *ncut, *degree, *xadj, *adjncy;
  float result;
  idxtype * adjwgt;

  ncut = idxsmalloc(npart, 0, "ComputeNCut: ncut");
  degree = idxsmalloc(npart, 0, "ComputeNCut: degree");
  if ( ncutVector == NULL )
  {
  	ncutVector=(float*)malloc(sizeof(float)*npart);
  }
  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  if (graph->adjwgt == NULL) {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
       		degree[cm] ++; 
        if (cm != where[adjncy[j]])
          ncut[cm] ++;
      }
    }
  }
  else {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
			degree[cm] += adjwgt[j];
        if (cm != where[adjncy[j]])
          ncut[cm] += adjwgt[j];
      }
    }
  }
  int empty = 0;
  result =0;
  for (i=0; i<npart; i++){
    if (degree[i] == 0)
      empty++;
    if (degree[i] >0)
	{
	  ncutVector[i] =ncut[i] *1.0/ degree[i]; 
      result += ncutVector[i];
	}
  }
  //printf("Empty clusters: %d\n", empty);
  free(ncut);
  free(degree);
  return result+empty;
}

/*************************************************************************
* This function computes the average conductance given the graph and a where vector
**************************************************************************/
float ComputeConductance(GraphType *graph, idxtype *where, int npart)
{
  int i, j, cm, nvtxs;
  idxtype *ncut, *degree, *xadj, *adjncy;
  float result;
  idxtype * adjwgt;

  ncut = idxsmalloc(npart, 0, "ComputeNCut: ncut");
  degree = idxsmalloc(npart, 0, "ComputeNCut: degree");
  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  idxtype totalGraphVol = 0;

  if (graph->adjwgt == NULL) {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
	       	degree[cm] ++; 
        if (cm != where[adjncy[j]])
          ncut[cm] ++;
      }
    }
  }
  else {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
			degree[cm] += adjwgt[j];
        if (cm != where[adjncy[j]])
          ncut[cm] += adjwgt[j];
      }
    }
  }

  for (i=0; i<npart; i++ )
  	totalGraphVol += degree[i];

  int empty = 0;
  result =0;
  for (i=0; i<npart; i++){
    if (degree[i] == 0)
      empty++;
    if (degree[i] >0)
	{
	  idxtype t = degree[i];
	  if ( degree[i] > totalGraphVol-degree[i] )
	  {
	  	t = totalGraphVol-degree[i];
		printf("Degree[i]:%d, totalGraphVol:%d\n", 
					degree[i], totalGraphVol );
	  }

      result +=  ncut[i] *1.0/t;
	  			
	}
  }
  //printf("Empty clusters: %d\n", empty);
  free(ncut);
  free(degree);
  return result+empty;
}



/*************************************************************************
* This function computes the normalized cut given the graph and a where vector
**************************************************************************/
float ComputeNCut(GraphType *graph, idxtype *where, int npart)
{
  int i, j, cm, nvtxs;
  idxtype *ncut, *degree, *xadj, *adjncy;
  float result;
  idxtype * adjwgt;

  ncut = idxsmalloc(npart, 0, "ComputeNCut: ncut");
  degree = idxsmalloc(npart, 0, "ComputeNCut: degree");
  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  if (graph->adjwgt == NULL) {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
	       	degree[cm] ++; 
        if (cm != where[adjncy[j]])
          ncut[cm] ++;
      }
    }
  }
  else {
    for (i=0; i<nvtxs; i++) {
      cm = where[i];
      for (j=xadj[i]; j<xadj[i+1]; j++){
	  	if ( adjncy[j] != i )
			degree[cm] += adjwgt[j];
        if (cm != where[adjncy[j]])
          ncut[cm] += adjwgt[j];
      }
    }
  }
  int empty = 0;
  result =0;
  for (i=0; i<npart; i++){
    if (degree[i] == 0)
      empty++;
    if (degree[i] >0)
      result +=  ncut[i] *1.0/ degree[i];
  }
  //printf("Empty clusters: %d\n", empty);
  free(ncut);
  free(degree);
 // return result+empty;
  return result;
}

int mapPartition(idxtype* part, idxtype nvtxs)
{
	int i,j,htcounter;
	idxtype* hashtable = idxmalloc(nvtxs,
							"hashtable for mapPartition");
	for(i=0;i<nvtxs;i++)
		hashtable[i]=-1;
	
	for(i=0,htcounter=0;i<nvtxs;i++)
	{
		if ( hashtable[part[i]] == -1 )
			hashtable[part[i]]=htcounter++;
	}

	for(i=0;i<nvtxs;i++)
		part[i]=hashtable[part[i]];

//	GKfree(&hashtable);
	free(hashtable);

	return htcounter;
}

/* Will assume that the the values are from 0 to numUnique-1 */
idxtype* histogram(idxtype* values, int n, int numUnique)
{
	idxtype* hist = idxmalloc(numUnique, "histogram:hist");
	int i;

	for ( i=0; i<numUnique; i++ )
		hist[i] = 0;

	for( i=0; i<n; i++ )
		hist[values[i]]++;
	
	return hist;
}

float stdDeviation(idxtype* values, int n)
{
	float mean = 0;
	double meanOfSquares = 0;

	int i;
	for ( i=0; i<n; i++ )
	{
		float t = values[i];
		mean += t;
		meanOfSquares += t*t; 
	}

	meanOfSquares /= n;
	mean /= n;
	if ( meanOfSquares < mean*mean )
	{
		printf("Yikes! meanOfSquares:%.4f, mean^2:%.4f",
			meanOfSquares, mean*mean);
		return 0;
	}

	return sqrt(meanOfSquares - mean*mean);
}

/*
int CheckBnd(GraphType *graph) 
{
  int i, j, nvtxs, nbnd;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    if (xadj[i+1]-xadj[i] == 0)
      nbnd++;   

    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) {
        nbnd++;
        ASSERT(bndptr[i] != -1);
        ASSERT(bndind[bndptr[i]] == i);
        break;
      }
    }
  }

  ASSERTP(nbnd == graph->nbnd, ("%d %d\n", nbnd, graph->nbnd));

  return 1;
}



int CheckBnd2(GraphType *graph) 
{
  int i, j, nvtxs, nbnd, id, ed;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    id = ed = 0;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      if (where[i] != where[adjncy[j]]) 
        ed += graph->adjwgt[j];
      else
        id += graph->adjwgt[j];
    }
    if (ed - id >= 0 && xadj[i] < xadj[i+1]) {
      nbnd++;
      ASSERTP(bndptr[i] != -1, ("%d %d %d\n", i, id, ed));
      ASSERT(bndind[bndptr[i]] == i);
    }
  }

  ASSERTP(nbnd == graph->nbnd, ("%d %d\n", nbnd, graph->nbnd));

  return 1;
}

int CheckNodeBnd(GraphType *graph, int onbnd) 
{
  int i, j, nvtxs, nbnd;
  idxtype *xadj, *adjncy, *where, *bndptr, *bndind;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;
  bndptr = graph->bndptr;
  bndind = graph->bndind;

  for (nbnd=0, i=0; i<nvtxs; i++) {
    if (where[i] == 2) 
      nbnd++;   
  }

  ASSERTP(nbnd == onbnd, ("%d %d\n", nbnd, onbnd));

  for (i=0; i<nvtxs; i++) {
    if (where[i] != 2) {
      ASSERTP(bndptr[i] == -1, ("%d %d\n", i, bndptr[i]));
    }
    else {
      ASSERTP(bndptr[i] != -1, ("%d %d\n", i, bndptr[i]));
    }
  }

  return 1;
}



int CheckRInfo(RInfoType *rinfo)
{
  int i, j;

  for (i=0; i<rinfo->ndegrees; i++) {
    for (j=i+1; j<rinfo->ndegrees; j++)
      ASSERTP(rinfo->edegrees[i].pid != rinfo->edegrees[j].pid, ("%d %d %d %d\n", i, j, rinfo->edegrees[i].pid, rinfo->edegrees[j].pid));
  }

  return 1;
}



int CheckNodePartitionParams(GraphType *graph)
{
  int i, j, k, l, nvtxs, me, other;
  idxtype *xadj, *adjncy, *adjwgt, *vwgt, *where;
  idxtype edegrees[2], pwgts[3];

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  where = graph->where;

  pwgts[0] = pwgts[1] = pwgts[2] = 0;
  for (i=0; i<nvtxs; i++) {
    me = where[i];
    pwgts[me] += vwgt[i];

    if (me == 2) { 
      edegrees[0] = edegrees[1] = 0;

      for (j=xadj[i]; j<xadj[i+1]; j++) {
        other = where[adjncy[j]];
        if (other != 2)
          edegrees[other] += vwgt[adjncy[j]];
      }
      if (edegrees[0] != graph->nrinfo[i].edegrees[0] || edegrees[1] != graph->nrinfo[i].edegrees[1]) {
        printf("Something wrong with edegrees: %d %d %d %d %d\n", i, edegrees[0], edegrees[1], graph->nrinfo[i].edegrees[0], graph->nrinfo[i].edegrees[1]);
        return 0;
      }
    }
  }

  if (pwgts[0] != graph->pwgts[0] || pwgts[1] != graph->pwgts[1] || pwgts[2] != graph->pwgts[2])
    printf("Something wrong with part-weights: %d %d %d %d %d %d\n", pwgts[0], pwgts[1], pwgts[2], graph->pwgts[0], graph->pwgts[1], graph->pwgts[2]);

  return 1;
}


int IsSeparable(GraphType *graph)
{
  int i, j, nvtxs, other;
  idxtype *xadj, *adjncy, *where;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  where = graph->where;

  for (i=0; i<nvtxs; i++) {
    if (where[i] == 2)
      continue;
    other = (where[i]+1)%2;
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      ASSERTP(where[adjncy[j]] != other, ("%d %d %d %d %d %d\n", i, where[i], adjncy[j], where[adjncy[j]], xadj[i+1]-xadj[i], xadj[adjncy[j]+1]-xadj[adjncy[j]]));
    }
  }

  return 1;
}
*/
int getLogBin(int a)
{
	if ( a <= 50 )
		return a;
	float la = log10(a);
	float la_trunc = round(la*100)/100.0;
	return (int)round(pow(10,la_trunc));
/*
	float flpart = la - floor(la);
	int base = (int)pow(10, floor(la));
	int ret = ;
	if ( flpart < 0.6 )
		ret = base + (base/10)*((int)round(flpart*100))
	else
		ret = base + base*((int)round(flpart*10));

	return ret;
*/	
}

idxtype* getWeightsHistogram(GraphType* graph, int* maxWeight, int
	logScale)
{
	int i;
	*maxWeight=0;
	int maxLogWeight;
	for ( i=0; i<graph->xadj[graph->nvtxs]; i++ )
	{
		if ( graph->adjwgt[i] > *maxWeight )
		{
			*maxWeight = graph->adjwgt[i];
			maxLogWeight = getLogBin(graph->adjwgt[i]);
		}
	}

	idxtype* hist;
	if ( logScale > 0 )
	{
		hist = idxsmalloc(maxLogWeight+1, 0,
					"getDegreeHistogram:hist");
	}
	else
	{
		hist = idxsmalloc(*maxWeight+1, 0,
							"getDegreeHistogram:hist");
	}

	for ( i=0; i<graph->xadj[graph->nvtxs]; i++ )
	{
		int l = graph->adjwgt[i];
		if ( logScale > 0 )
		{
			l = getLogBin(l);
		}
		hist[l]++;
	}
	
	return hist;
}

idxtype* getDegreeHistogram(GraphType* graph, int* maxDegree, int
	logScale)
{
	int i;
	*maxDegree=0;
	int maxLogDegree;
	for ( i=0; i<graph->nvtxs; i++ )
	{
		int k;
		if ( (k=(graph->xadj[i+1] - graph->xadj[i])) > *maxDegree )
		{
			*maxDegree = k;
			maxLogDegree = getLogBin(k);
		}
	}

	idxtype* hist;
	if ( logScale > 0 )
	{
		hist = idxsmalloc(maxLogDegree+1, 0,
					"getDegreeHistogram:hist");
	}
	else
	{
		hist = idxsmalloc(*maxDegree+1, 0,
							"getDegreeHistogram:hist");
	}

	for ( i=0; i<graph->nvtxs; i++ )
	{
		int l = graph->xadj[i+1]-graph->xadj[i];
		if ( logScale > 0 )
		{
			l = getLogBin(l);
		}
		hist[l]++;
	}
	
	return hist;
}

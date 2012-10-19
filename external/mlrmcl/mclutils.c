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



#include <metis.h>

void dump_graph(GraphType *graph)
{
	int i,j;
	int nvtxs = graph->nvtxs;
	int wgtflag = (graph->adjwgt==NULL)?0:1;
	FILE *fp=stderr;
	
	for (i=0; i<nvtxs; i++)
	{
		fprintf(fp, "%d", i);
		for (j = graph->xadj[i];j < graph->xadj[i+1];j++)
		{
			fprintf(fp, " %d", graph->adjncy[j]);
			if ( wgtflag )
				fprintf(fp, ":%d", graph->adjwgt[j]);
		}
		fprintf(fp, "\n");
	}

	fprintf(fp, "\n");
/*	if (graph->rmap1 != NULL )
	{
		for( i=0; i<nvtxs; i++ )
		{
			fprintf(fp,"%d <- %d", i, graph->rmap1[i]);
			if (graph->rmap2[i] > -1)
				fprintf(fp, ",%d", graph->rmap2[i]);
			fprintf(fp, "\n");
		}
	}
*/
}


void testCoarsening(int *nvtxs, idxtype *xadj, idxtype *adjncy,
		idxtype *vwgt, idxtype *adjwgt, 
		int *wgtflag, int ct)
{
	GraphType graph, *cgraph;
	CtrlType ctrl;

	my_SetUpGraph(&graph, *nvtxs, xadj, adjncy, vwgt, adjwgt,
	*wgtflag, 1);
	/* The last argument indicates we are setting up the original
	 * graph */

	ctrl.CoarsenTo=ct;
	ctrl.CType=MATCH_SHEMN;

	my_AllocateWorkSpace(&ctrl,&graph);

	cgraph = Coarsen2Way(&ctrl,&graph);	

	do{
		dump_graph(cgraph);
		cgraph = cgraph->finer;
	}while ( cgraph != NULL );

}

void ncutifyWeights(Matrix* M, int hasWeightsAlready, int ncutify)
{
	int i,j;
	if ( hasWeightsAlready )
	{
		for(i=0;i<M->nvtxs;i++)
		{
			if ( ncutify==4)
				M->adjwgtsum[i]=(wgttype)1.0/sqrt(M->adjwgtsum[i]);
			else
				M->adjwgtsum[i] = 1.0/M->adjwgtsum[i];
			if (ncutify==2)
				M->adjwgtsum[i] *= M->adjwgtsum[i];
		}
	}
	else
	{
		for(i=0;i<M->nvtxs;i++)
		{
			if ( ncutify==4)
				M->adjwgtsum[i]=(wgttype)1.0/sqrt((wgttype)(M->xadj[i+1]-M->xadj[i]));
			else
				M->adjwgtsum[i] = 1.0/(M->xadj[i+1]-M->xadj[i]);
			if (ncutify==2)
				M->adjwgtsum[i]=M->adjwgtsum[i]*M->adjwgtsum[i];
		}
	}

	if ( hasWeightsAlready )
	{
		for(i=0;i<M->nvtxs;i++)
		{
			for(j=M->xadj[i];j<M->xadj[i+1];j++)
			{
				if (ncutify==4)
					M->adjwgt[j]*=M->adjwgtsum[i]*M->adjwgtsum[M->adjncy[j]];
				else
					M->adjwgt[j] *=
						(M->adjwgtsum[i] + M->adjwgtsum[M->adjncy[j]]);
				if ( ncutify==3)
					M->adjwgt[j] *= M->adjwgt[j];
			}
		}
	}
	else
	{
		for(i=0;i<M->nvtxs;i++)
		{
			for(j=M->xadj[i];j<M->xadj[i+1];j++)
			{
				if (ncutify==4)
					M->adjwgt[j]=M->adjwgtsum[i]*M->adjwgtsum[M->adjncy[j]];
				else
					M->adjwgt[j] =
					(M->adjwgtsum[i] + M->adjwgtsum[M->adjncy[j]]);
				if ( ncutify==3)
					M->adjwgt[j] *= M->adjwgt[j];
			}
		}
	}
//	printf("sample wgt after ncutify:%1.3f\n",M->adjwgt[0]);
}

void normalizeColumns(Matrix* M, int hasWeightsAlready, int
hasSumsAlready)
{
	int i,j;
	wgttype sum;
	for(i=0; i<M->nvtxs; i++)
	{
		sum=0;
		if ( hasWeightsAlready )
		{
			if ( !hasSumsAlready )
			{
				for(j=M->xadj[i];j<M->xadj[i+1];j++)
				{
					sum+=M->adjwgt[j];
				}
			}
			else
				sum=M->adjwgtsum[i];
		}
		else
			sum=(wgttype)(M->xadj[i+1]-M->xadj[i]);

		if ( hasWeightsAlready )
		{
			for(j=M->xadj[i];j<M->xadj[i+1];j++)
				M->adjwgt[j] /= sum;
		}
		else
		{
			for(j=M->xadj[i];j<M->xadj[i+1];j++)
				M->adjwgt[j] = 1.0/sum;
		}


		if ( M->adjwgtsum != NULL )
			M->adjwgtsum[i]=1.0;
	}
//	printf("sample weight after normalize:%1.3f\n",M->adjwgt[0]);
}

void sortAdjLists(int nvtxs, idxtype* xadj, idxtype* adjncy, wgttype* adjwgt)
{
	int i,j;
	for(i=0;i<nvtxs;i++)
	{
		ParallelQSort(adjncy,adjwgt,xadj[i],xadj[i+1]-1);
	}
}

/* NOTE: I assume there are no self loops already. */
Matrix* addSelfLoops(Matrix* M)
{
	Matrix* ret=allocMatrix(M->nvtxs, M->nnz+M->nvtxs, 0, 0, 0);
	int i=0;

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int j,jj;
		/* NOTE: I assume the adjacency lists are sorted in increasing
		 * order! */
		int addedLoop=0;
		for ( j=M->xadj[i], jj=ret->xadj[i]; j<M->xadj[i+1]; j++)
		{
			if ( M->adjncy[j] > i && !addedLoop )
			{
				ret->adjncy[jj] = i;
				ret->adjwgt[jj++] = 1;
				addedLoop=1;
			}
			ret->adjncy[jj] = M->adjncy[j];
			ret->adjwgt[jj++] = M->adjwgt[j];
		}
		if ( !addedLoop )
		{
			ret->adjncy[jj] = i;
			ret->adjwgt[jj++] = 1;
			addedLoop=1;
		}
		ret->xadj[i+1] = jj;
	}

	return ret;
}

void removeSelfLoops(GraphType* g)
{
	int i,j,k;

	j=0;
	for ( i=0; i<g->nvtxs; i++)
	{
		k = g->xadj[i];
		g->xadj[i] = j;
		for ( ; k<g->xadj[i+1]; k++ )
		{
			if ( g->adjncy[k] == i )
				continue;
			g->adjncy[j] = g->adjncy[k];
			if ( g->adjwgt != NULL )
				g->adjwgt[j] = g->adjwgt[k];
			j++;
		}
	}
	g->xadj[g->nvtxs] = j;
	g->nedges = j;
}

Matrix* removeSelfLoops(Matrix* M, int normalizeColumns)
{
	int nvtxs=M->nvtxs,nnz=M->nnz,i, ret_counter;
	Matrix* ret=allocMatrix(nvtxs, nnz, 0, 0, 0);
	wgttype sum;

	ret->xadj[0]=0;
	ret_counter=0;
	for ( i=0; i<nvtxs; i++ )
	{
		int j;
		sum=0;
		for( j=M->xadj[i]; j<M->xadj[i+1]; j++ )
		{
			if ( M->adjncy[j] == i )
				continue;
			ret->adjncy[ret_counter]=M->adjncy[j];
			ret->adjwgt[ret_counter++]=M->adjwgt[j];
			sum+=M->adjwgt[j];
		}
		ret->xadj[i+1]=ret_counter;

		if ( normalizeColumns > 0 && sum > 0 )
		{
			for( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
				ret->adjwgt[j]/=sum;
		}
	}

	return ret;
}

Matrix* setupCanonicalMatrix(int nvtxs, int nedges, idxtype* xadj,
	idxtype* adjncy, idxtype* adjwgt, int ncutify)
{
	int i,j;
	Matrix* ret;
	if ( ncutify )
		ret=allocMatrix(nvtxs,nedges,1,0,0);
	else
		ret=allocMatrix(nvtxs,nedges,0,0,0);
		
	idxcopy(nvtxs+1, xadj, ret->xadj); 
	idxcopy(nedges, adjncy, ret->adjncy);
	if ( adjwgt != NULL )
	{
		if ( ncutify )
		{
			for(i=0;i<ret->nvtxs;i++)
			{
				ret->adjwgtsum[i]=0;
				for(j=ret->xadj[i];j<ret->xadj[i+1];j++)
				{
					ret->adjwgt[j]=(wgttype)adjwgt[j];
					ret->adjwgtsum[i]+=ret->adjwgt[j];
				}
			}
			ncutifyWeights(ret,1,ncutify);
		}
		else
		{
			for(i=0;i<nedges;i++)
				ret->adjwgt[i]=(wgttype)adjwgt[i];
		}
		normalizeColumns(ret,1,0);
	}
	else
	{
		if ( ncutify )
			ncutifyWeights(ret,0,ncutify);
		normalizeColumns(ret,0,0);
	}

	// sort each column in ascending order. This is necessary for
	// getDprAdjMatrix. 
	for(i=0;i<nvtxs;i++)
	{
		ParallelQSort(ret->adjncy,ret->adjwgt,ret->xadj[i],ret->xadj[i+1]-1);
	}

	return ret;
}

void printRow(Matrix *M, int row)
{
	int i,j;
	wgttype sum=0;
	printf("Row of %d:", row);
	for(i=0;i<M->nvtxs;i++)
	{
		for(j=M->xadj[i];j<M->xadj[i+1];j++)
		{
			if ( M->adjncy[j] == row )
			{
				printf(" %d:%1.3f",i,M->adjwgt[j]);
				sum += M->adjwgt[j];
			}
		}
	}
	printf("\n");
	printf("Sum of row:%f\n", sum);
}

int numWithoutAttractors(Matrix* M)
{
	int i,numWithoutAttrs=0;
	for(i=0;i<M->nvtxs;i++)
	{
		if (M->attractors[i] == -1)
			numWithoutAttrs++;
	}
	return numWithoutAttrs;
}

void getAttractorsForAll(Matrix* M)
{
	/* A node is assigned an attractor in pruneAndNormalize only
	 * if there's a node to which it has a flow more than 0.5, 
	 * otherwise the attractor is -1. This is probably how it
	 * should be because the attractors computed in
	 * pruneAndNormalize are used to check for convergence. Now,
	 * once we've decided that we've converged, we need to assign
	 * an attractor to a node, even if the max. flow in that node
	 * is less than 0.5. We simply assign the node with the max
	 * flow as the attractor to that node. */

	 if ( M->attractors == NULL )
	 {
	 	M->attractors=idxmalloc(M->nvtxs,
						"getAttractorsForAll:M->attractors");
	 }

	 int i,j,attr,numWithoutAttrs=0, numTornNodes=0;
	 wgttype fl_maxwgt;
	 int avgTornNodeDegree=0;
	 int int_maxwgt;
	 for(i=0;i<M->nvtxs;i++)
	 {
	 	if (M->attractors[i] != -1 )
			continue;
		fl_maxwgt = 0;
		int_maxwgt = 0;
		attr=-1;
		for(j=M->xadj[i];j<M->xadj[i+1];j++)
		{
			if ( M->adjwgt[j] > fl_maxwgt )
			{
				attr=M->adjncy[j];
				fl_maxwgt = M->adjwgt[j];
				int_maxwgt = M->adjwgt[j];
			}
		}
		M->attractors[i]=attr;
		if ( fl_maxwgt < 0.95 )
		{
			numWithoutAttrs++;
/*			if ( i < 50 )
			{
				printf("maxwgt of node %d:%f\n", i+1, fl_maxwgt);
				printf("degree: %d\n", M->xadj[i+1]-M->xadj[i]);
			}
*/		}
		float diff = fabs(fl_maxwgt -
		1.0/(M->xadj[i+1]-M->xadj[i]));
		if (  diff < 0.01 )
		{
//			printf("torn node: %d\n", i+1);
			if ( i < 50 )
			{
//				printf("diff:%f\n",diff ); 
			}
			numTornNodes++;
			avgTornNodeDegree += M->xadj[i+1]-M->xadj[i];
		}

	 }
//	 printf("Number of nodes without attractors:%d\n",numWithoutAttrs);
//	 printf("Number of torn nodes:%d\n", numTornNodes);
//	 printf("Avg. torn node degree:%f\n",
//			 avgTornNodeDegree/((float)numTornNodes) );
}

int isConverged(idxtype* oldAttr, idxtype* newAttr, int n)
{
	/* For now, I am just comparing attractors. But the issue of
	 * convergence seems to need some serious thought. What to do
	 * with nodes which have a foot each in multiple clusters?
	 * They may not have an attractor (a row with >.5 weight), and
	 * even if they do, that may keep fluctuating. */
	int i;
	if ( oldAttr == NULL || newAttr == NULL )
		return 0;

	for( i=0; i<n; i++ )
	{
/*		printf("%d: oldAttr:%d,newAttr:%d\n",
					i,oldAttr[i],newAttr[i]);
*/		if( oldAttr[i]<0 || newAttr[i]<0 || newAttr[i] != oldAttr[i] )
			return 0;
	}

	return 1;
}

int isConverged2(idxtype* oldAttr, idxtype* newAttr, int n)
{
	/* For now, I am just comparing attractors. But the issue of
	 * convergence seems to need some serious thought. What to do
	 * with nodes which have a foot each in multiple clusters?
	 * They may not have an attractor (a row with >.5 weight), and
	 * even if they do, that may keep fluctuating. */
	int i,diff=0;
	if ( oldAttr == NULL || newAttr == NULL )
		return n;

	for( i=0; i<n; i++ )
	{
/*		printf("%d: oldAttr:%d,newAttr:%d\n",
					i,oldAttr[i],newAttr[i]);
*/		if(/*oldAttr[i]<0 || newAttr[i]<0 ||*/newAttr[i] != oldAttr[i] )
			diff++;		
	}

	return diff;
}

int noOfUnique(idxtype* attrs,int n, idxtype* hashtable)
{
	// will assume that the hashtable is clear, i.e. that all
	// entries are -1. 
	int htcounter=0,i=0;
	for(i=0;i<n;i++)
	{
	/*	if ( attrs[i] < 0 || attrs[i] >= n )
		{
			printf("attrs[%d]:%d\n",i,attrs[i]);
			abort();
		} */
		if ( attrs[i]>=0 && hashtable[attrs[i]] < 0)
		{
			hashtable[attrs[i]]=1;
			htcounter++;
		}
	}
	for(i=0;i<n;i++)
	{
		hashtable[i]=-1;
	}
	return htcounter;
}

void freeMatrix(Matrix* a)
{
	if ( a == NULL )
		return;

	if ( a->adjncy != NULL )
	{
		free(a->adjncy);
		free(a->xadj);
		free(a->adjwgt);
	}
	if ( a->adjwgtsum != NULL )
		free(a->adjwgtsum);
	if ( a->maxwgt != NULL )
		free(a->maxwgt);
	if ( a->attractors != NULL )
		free(a->attractors);
	if ( a->rmap != NULL )
		free(a->rmap);

	free(a);
	a=NULL;
}

Matrix* allocMatrix(int nvtxs, int nedges, int allocSum, int
						allocMax, int allocAttractor)
{
	Matrix* a;
	a=(Matrix*)malloc(sizeof(Matrix));
	a->nvtxs=nvtxs;
	a->nnz=nedges;
	a->xadj = (idxtype*) GKmalloc( sizeof(idxtype)*(nvtxs+1) , 
					"allocMatrix:xadj" );
	a->adjncy=(idxtype*) GKmalloc(sizeof(idxtype)*nedges , 
						"allocMatrix:adjncy" );
	a->adjwgt=(wgttype*)GKmalloc( sizeof(wgttype)*nedges , 
						"allocMatrix:adjwgt" );
	if ( allocSum )
		a->adjwgtsum=(wgttype*)malloc(sizeof(wgttype)*nvtxs);
	else
		a->adjwgtsum=NULL;
		
 	if ( allocSum )
		a->maxwgt=(wgttype*)malloc(sizeof(wgttype)*nvtxs);
	else
		a->maxwgt=NULL;

	if ( allocAttractor )
	{
		a->attractors=(idxtype*)malloc(sizeof(idxtype)*nvtxs);
//		printf("Allocating attractors\n");
	}
	else
		a->attractors=NULL;

	a->rmap=NULL;
	a->currentSize = nedges;
	return a;
}

void dumpMatrix(Matrix* a)
{
	FILE* fp=stderr;
	int i,j;
	for( i=0; i<a->nvtxs;i++)
	{
		fprintf(fp,"%d",i);
		for(j=a->xadj[i];j<a->xadj[i+1];j++)
			fprintf(fp," %d:%1.4f", a->adjncy[j],a->adjwgt[j]);
		fprintf(fp,"\n");
	}
}

void getPermutedGraph(idxtype* perm, idxtype* revPerm, int nvtxs, 
		int nedges, idxtype* xadj, idxtype* adjncy, idxtype* adjwgt, 
		idxtype** p_xadj, idxtype** p_adjncy, idxtype** p_adjwgt)
{
	*p_xadj = idxmalloc(nvtxs+1, "getPermutedGraph:p_xadj");
	*p_adjncy = idxmalloc(nedges, "getPermutedGraph:p_adjncy");
	if ( adjwgt != NULL )
		*p_adjwgt = idxmalloc(nedges, "getPermutedGraph:p_adjwgt");
	else
		*p_adjwgt = NULL;

	int i;
	(*p_xadj)[0]=0;
	for ( i=0; i<nvtxs; i++ )
	{
		int orgI = revPerm[i];
		int j;
		(*p_xadj)[i+1] = (*p_xadj)[i] + ( xadj[orgI+1] -
								xadj[orgI] );
		for ( j=(*p_xadj)[i]; j<(*p_xadj)[i+1]; j++ )
		{
			int orgJ = xadj[orgI] + j - (*p_xadj)[i];
			(*p_adjncy)[j] = perm[adjncy[orgJ]];
			if ( adjwgt != NULL )
				(*p_adjwgt)[j] = adjwgt[orgJ];
		}
		
		if ( adjwgt != NULL )
		{
			ParallelQSortInts( (*p_adjncy), (*p_adjwgt), (*p_xadj)[i],
			(*p_xadj)[i+1]-1 );						
		}
		else
		{
			iidxsort( (*p_xadj)[i+1]-(*p_xadj)[i],
						(*p_adjncy)+(*p_xadj)[i] );	
		}
	}

	return;
}

// this actually undoes the permutation
Matrix* permuteRowsAndColumns(Matrix* M, idxtype* perm)
{
	Matrix *ret = allocMatrix(M->nvtxs, M->nnz, 0, 0, 0);
	int i;
	idxtype* revPerm = idxmalloc(M->nvtxs,
						"permuteRowsAndColumns:revPerm");

	for ( i=0; i<M->nvtxs; i++ )
	{
		revPerm[perm[i]] = i;
	}

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int orgI = revPerm[i];
		int j;
		ret->xadj[i+1] = ret->xadj[i] + ( M->xadj[orgI+1] -
								M->xadj[orgI] );
		for ( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
		{
			int orgJ = M->xadj[orgI] + j - ret->xadj[i];
			ret->adjncy[j] = perm[M->adjncy[orgJ]];
			ret->adjwgt[j] = M->adjwgt[orgJ];
		}
								
		ParallelQSort( ret->adjncy, ret->adjwgt, ret->xadj[i],
		ret->xadj[i+1]-1 );						
	}

	ret->nnz = M->nnz;

	free ( revPerm );

	return ret;
}

// this actually undoes the permutation
Matrix* permuteRowsAndColumns(Matrix *M, idxtype* rowPerm, idxtype* colPerm)
{
	Matrix *ret = allocMatrix(M->nvtxs, M->nnz, 0, 0, 0);
	int i;
	idxtype* revRowPerm = idxmalloc(M->nvtxs,
						"permuteRowsAndColumns:revRowPerm");
//	idxtype* revColPerm = idxmalloc(M->nvtxs, 
//						"permuteRowsAndColumns:revColPerm");

	for ( i=0; i<M->nvtxs; i++ )
	{
		revRowPerm[rowPerm[i]] = i;
//		revColPerm[colPerm[i]] = i;
	}

	ret->xadj[0]=0;
	for ( i=0; i<M->nvtxs; i++ )
	{
		int orgI = revRowPerm[i];
		int j;
		ret->xadj[i+1] = ret->xadj[i] + ( M->xadj[orgI+1] -
								M->xadj[orgI] );
		for ( j=ret->xadj[i]; j<ret->xadj[i+1]; j++ )
		{
			int orgJ = M->xadj[orgI] + j - ret->xadj[i];
			ret->adjncy[j] = colPerm[M->adjncy[orgJ]];
			ret->adjwgt[j] = M->adjwgt[orgJ];
		}
								
		ParallelQSort( ret->adjncy, ret->adjwgt, ret->xadj[i],
		ret->xadj[i+1]-1 );						
	}

	ret->nnz = M->nnz;

//	GKfree ( (void**)&revRowPerm, (void**)&revColPerm, LTERM );
	GKfree ( (void**)&revRowPerm, LTERM );

	return ret;
	
}

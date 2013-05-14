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


/** This file will have routines for obtaining a subgraph of the
 * present graph induced by a subset of nodes of the original
 * graph. For now, we won't deal with issues of the
 * induced subgraph being disconnected; the input subset is
 * hopefully nice enough that the induced subgraph won't get
 * disconneted.
 */

#include <metis.h>

void getSubgraph(GraphType *graph, idxtype* newIds, int new_nvtxs, int
wgtflag, GraphType **new_graph)
{
	int i=0, j=0, new_nedges=0;
	idxtype *new_xadj, *new_adjncy, *new_adjwgt, *new_vwgt;

	for ( i=0, new_nedges=0; i<graph->nvtxs; i++ )
	{
		if ( newIds[i] > -1 )
			new_nedges += graph->xadj[i+1]-graph->xadj[i];
	}

	new_xadj = idxmalloc(new_nvtxs + 1,
								"getSubgraph:new_xadj");
	new_adjncy = idxmalloc(new_nedges,"getSubgraph: new_adjncy");
	
	if ( graph->adjwgt != NULL )
	{
		new_adjwgt = idxmalloc(new_nedges, "getSubgraph: new_adjwgt");
	}

	if ( graph->vwgt != NULL )
	{
		new_vwgt = idxmalloc(new_nvtxs, "getSubgraph: new_vwgt");
		for( i=0; i<graph->nvtxs; i++ )
		{
			if ( newIds[i] > -1 )
				new_vwgt[newIds[i]] = graph->vwgt[i];
		}
	}

	new_xadj[0] = 0;
	int edgeCounter=0;
	for ( i=0; i<graph->nvtxs; i++ )
	{
		if ( newIds[i] < 0 )
			continue;
		for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
		{
			if ( newIds[graph->adjncy[j]] < 0 )
				continue;
			new_adjncy[edgeCounter] = newIds[graph->adjncy[j]];
			if ( graph->adjwgt != NULL )
			{
				new_adjwgt[edgeCounter] = graph->adjwgt[j];
			}
			edgeCounter++;
		}
		new_xadj[newIds[i]+1] = edgeCounter;
	}

//	free(graph->xadj);
//	free(graph->adjncy);

	*new_graph=(GraphType*)malloc(sizeof(GraphType));
	my_SetUpGraph(*new_graph, new_nvtxs, new_xadj, new_adjncy,
	new_vwgt, new_adjwgt, wgtflag, 1);

//	FreeGraph(graph);
//	free(graph->gdata);
//	free(graph);
//	graph = new_graph;
	
//	return newIds;
}	

/** 
 * Removes nodes with degree greater than threshold from the
 * graph.
 * Returns the mapping of old node ids to new node ids.
 */
idxtype* removeHubs(GraphType *graph, int threshold, int wgtflag,
GraphType **new_graph, int removeSingletons)
{
	int i=0, j=0, numHubs=0, numSingletons;
	idxtype *newIds;
	
	for ( i=0; i<graph->nvtxs; i++ )
	{
		if ( graph->xadj[i+1]-graph->xadj[i] > threshold )
			numHubs++;
		else if ( removeSingletons > 0 )
		{
			if ( graph->xadj[i+1]==graph->xadj[i] )
				numSingletons++;
			else if ( graph->xadj[i+1]==graph->xadj[i]+1 
			   && graph->adjncy[graph->xadj[i]] == i )
			   	numSingletons++;
		}
	}

	printf("Num hubs to remove:%d\n", numHubs);
	if ( removeSingletons > 0 )
		printf("Num singletons to remove:%d\n", numSingletons);

	//listOfHubs = idxmalloc(numHubs, "removeHubs:listOfHubs");
	newIds = idxmalloc(graph->nvtxs, "removeHubs:newIds");

	int newIdCounter=0;
	for ( i=0; i<graph->nvtxs; i++ )
	{
		newIds[i]=newIdCounter++;
		if ( graph->xadj[i+1]-graph->xadj[i] > threshold )
		{
			newIds[i]=-1;
			newIdCounter--;
	//		listOfHubs[hubCounter++]=i;
//			totalHubDegree += graph->xadj[i+1]-graph->xadj[i];
		}
		else if ( removeSingletons > 0 )
		{
			if ( graph->xadj[i+1]==graph->xadj[i] 
			 || ( graph->xadj[i+1]==graph->xadj[i]+1 
			   && graph->adjncy[graph->xadj[i]] == i ) )
			{
				newIds[i]=-1;
				newIdCounter--;
			}
		}
	}

	getSubgraph(graph, newIds, newIdCounter, wgtflag, new_graph);

	return newIds;
}

GraphType* getCutGraph(GraphType* graph, idxtype* indices, 
						int wgtflag)
{
	int i,j,k;
	GraphType* cutGraph = (GraphType*)malloc(sizeof(GraphType));

	int nedges = ComputeNumEdgesCut(graph, indices);
	idxtype *adjncy, *adjwgt, *xadj;
	cutGraph->adjncy = idxmalloc( nedges, "getCutGraph:cutGraph->adjncy");
	cutGraph->xadj = idxmalloc( graph->nvtxs, "getCutGraph:cutGraph->adjwgt");
	if ( graph->adjwgt != NULL )
		cutGraph->adjwgt = idxmalloc( nedges, "getCutGraph:cutGraph->adjwgt");
	
	cutGraph->xadj[0] = 0;
	for ( i=0, k=0; i<graph->nvtxs; i++ )
	{
		for ( j=graph->xadj[i]; j<graph->xadj[i+1]; j++ )
		{
			int b = graph->adjncy[j];
			if ( indices[i] != indices[b] )
			{
				cutGraph->adjncy[k] = b;
				if ( graph->adjwgt != NULL )	
					cutGraph->adjwgt[k] = graph->adjwgt[j];
				k++;
			}
		}
		cutGraph->xadj[i+1] = k;
	}

	return cutGraph;
}

void globallySampleEdges(int nvtxs, int nedges, idxtype* xadj, idxtype*
	adjncy, idxtype* adjwgt, idxtype** new_xadj, idxtype**
	new_adjncy, float gsRatio)
{
	int numEdgesToSample = (int) round(gsRatio*nedges);
	printf("numEdgesToSample:%d\n", numEdgesToSample);

	int *edgeIds = idxmalloc(numEdgesToSample, 
						"globallySampleEdges:edgeIds");

	int *edgeSampled = idxmalloc(nedges, 
					"globallySampleEdges:edgesSampled");

	int i;
	for ( i=0; i<nedges; i++ )
		edgeSampled[i] = 0;

	for ( i=0; i<numEdgesToSample; i++ )
	{
		int rnd;
		do {
			rnd = RandomInRangeFast(nedges);
		} while ( edgeSampled[rnd] > 0 );

		edgeIds[i] = rnd;
		edgeSampled[rnd] = 1;
	}

	iidxsort(numEdgesToSample, edgeIds);

	(*new_xadj) = idxmalloc(nvtxs+1, "globallySampleEdges:new_xadj");
	(*new_adjncy) = idxmalloc(numEdgesToSample,
					"globallySampleEdges:new_adjncy");

	(*new_xadj)[0]=0;
	int j;
	for ( i=0, j=0; i<numEdgesToSample; i++ )
	{
		while ( edgeIds[i] >= xadj[j+1] )
		{
			(*new_xadj)[j+1]=i;
			j++;
		}

		(*new_adjncy)[i] = adjncy[edgeIds[i]]; 
	}

	while ( j < nvtxs )
	{
		(*new_xadj)[j+1] = i;
		j++;
	}

	free(edgeIds);
	free(edgeSampled);
}


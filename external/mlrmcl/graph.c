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
 * $Id: graph.c,v 1.3 2010/03/12 00:25:53 venu Exp $
 *
 */

#include <metis.h>


void my_SetUpGraph(GraphType *graph, int nvtxs, idxtype* xadj,
idxtype* adjncy, idxtype* vwgt, idxtype* adjwgt, int
wgtflag, int isOrgGraph )
{
	int i, sum, j, k, gsize;
	InitGraph(graph);
	
	graph->nvtxs = nvtxs;
    graph->nedges = xadj[nvtxs];
    graph->ncon = 0;
    graph->xadj = xadj;
    graph->adjncy = adjncy;
	graph->numDescendants = NULL;
	graph->vwgt = NULL;
	
    gsize = 0; 
    if ((wgtflag&2) == 0)
      gsize += nvtxs;
    if ((wgtflag&1) == 0)
      gsize += graph->nedges;

    gsize += 2*nvtxs;
	/* need extra space for storing the rmaps, but we don't
	 * need rmaps for the original graph. */
	if ( isOrgGraph != 1 )
    	gsize += 2*nvtxs; 

    graph->gdata = idxmalloc(gsize, "SetUpGraph: gdata");

    /* Create the vertex/edge weight vectors if they are not supplied */
    gsize = 0;
    if ((wgtflag&2) == 0) {
      vwgt = graph->vwgt = idxset(nvtxs, 1, graph->gdata);
      gsize += nvtxs;
    }
    else
      graph->vwgt = vwgt;

    if ((wgtflag&1) == 0) {
      adjwgt = graph->adjwgt = idxset(graph->nedges, 1, graph->gdata+gsize);
      gsize += graph->nedges;
    }
    else
      graph->adjwgt = adjwgt;


    /* Compute the initial values of the adjwgtsum */
    graph->adjwgtsum = graph->gdata + gsize;
    gsize += nvtxs;

	ComputeAdjWgtSums(graph);
/*    for (i=0; i<nvtxs; i++) {
      sum = 0;
      for (j=xadj[i]; j<xadj[i+1]; j++)
        sum += adjwgt[j];
      graph->adjwgtsum[i] = sum;
    }
*/
    graph->cmap = graph->gdata + gsize;
    gsize += nvtxs;

	if ( isOrgGraph != 1 )
	{
		graph->rmap1 = graph->gdata + gsize;
		gsize += nvtxs;
//		graph->rmap2 = graph->gdata + gsize;
//		gsize += nvtxs;
	}
	else
	{
		graph->rmap1=graph->rmap2=NULL;
	}

	graph->isOrgGraph=isOrgGraph;
	graph->coarser = NULL;
	graph->finer = NULL;

}


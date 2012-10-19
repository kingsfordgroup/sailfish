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
 * ccgraph.c
 *
 * This file contains the functions that create the coarse graph
 *
 * Started 8/11/97
 * George
 *
 * $Id: ccgraph.c,v 1.6 2010/01/24 06:23:16 venu Exp $
 *
 */

#include <metis.h>


/* Venu: my version of createcoarsegraph */

void my_CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
	CreateCoarseGraph(ctrl,graph,cnvtxs,match,perm);
}
/*  int i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask; 
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *rmap1, *rmap2, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  mask = HTLENGTH;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) { 
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;
  rmap1 = graph->rmap1;
  rmap2 = graph->rmap2;
*/
  /* Initialize the coarser graph */
/*  cgraph = my_SetUpCoarseGraph(graph, cnvtxs);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;

  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0]=0;
  for(i=0,ii=0; i<cnvtxs;i++)
  {
  	l = rmap1[i];
	m = rmap2[i];

	for( jj=xadj[l]; jj<xadj[l+1]; jj++)
		{
			cadjncy[ii]=cmap[adjncy[jj]];
			if ( adjwgt = NULL )
				cadjwgt[ii]=1;
			else
				cadjwgt[ii]=adjwgt[jj];
		}
	if ( m > -1 )
	{
		for( jj=xadj[m];jj<xadj[m+1];jj+)
		{
			
		}
	}
  }
  cgraph->nedges = cnedges;

*/
  /* Ditch readjustmemory. In any case the code in this function
   * seems to be inside an if condition that will never be
   * satisfied. */
/*  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}
*/

void CreateCoarseGraph_PowerLaw(CtrlType *ctrl, GraphType *graph,
		int cnvtxs, idxtype *match, idxtype *cpointers, idxtype
		*cvwgt)
{
	int i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask, dovsize;
	idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
	idxtype *cmap, *htable;
	idxtype *cxadj, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
	float *nvwgt, *cnvwgt;
	GraphType *cgraph;

	dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

	mask = HTLENGTH;

	// Ditch this optimization. 
/*	if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) 
	{ 
		CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
		return;
	}
*/
	IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

	nvtxs = graph->nvtxs;
	ncon = graph->ncon;
	xadj = graph->xadj;
	vwgt = graph->vwgt;
	vsize = graph->vsize;
	nvwgt = graph->nvwgt;
	adjncy = graph->adjncy;
	adjwgt = graph->adjwgt;
	adjwgtsum = graph->adjwgtsum;
	cmap = graph->cmap;

	/* Initialize the coarser graph */
	/* cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
	*/
	/* Venu: my addition */
	cgraph = my_SetUpCoarseGraph(graph, cnvtxs);
	cxadj = cgraph->xadj;
	idxcopy(cnvtxs, cvwgt, cgraph->vwgt);
	cvsize = cgraph->vsize;
	cnvwgt = cgraph->nvwgt;
	cadjwgtsum = cgraph->adjwgtsum;
	cadjncy = cgraph->adjncy;
	cadjwgt = cgraph->adjwgt;

	iend = xadj[nvtxs];
	auxadj = ctrl->wspace.auxcore; 
	memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 

	for (i=0; i<iend; i++)
		auxadj[i] = cmap[auxadj[i]];

	htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

	float maxvwgt=0, avgvwgt=0;
	for ( i=0; i<cnvtxs; i++ )
	{
		avgvwgt += cgraph->vwgt[i];
		if ( cgraph->vwgt[i] > maxvwgt )
			maxvwgt = cvwgt[i];
	}
	avgvwgt /= cnvtxs;
	printf("Max. cvwgt: %.1f, Avg. cvwgt:%.2f\n",
			 maxvwgt, avgvwgt);
	fflush(stdout);


	cxadj[0] = cnedges = 0;
	for ( i=0; i<cnvtxs; i++ )
	{
		for ( j=cpointers[i]; j != UNMATCHED; j=match[j])
		{
			if ( xadj[j+1] == xadj[j] )
			{
				printf("Zero neighbors for j:%d\n", j);
				abort();
			}
			for ( jj=xadj[j]; jj<xadj[j+1]; jj++ )
			{
				int k=auxadj[jj];
				int kk = k&mask;

				if ( (m=htable[kk]) == -1 )
				{
					cadjncy[cnedges] = k;
					cadjwgt[cnedges] = adjwgt[jj];
					htable[kk] = cnedges++;
				}
				else if ( cadjncy[m] == k )
				{
					cadjwgt[m] += adjwgt[j];
				}
				else
				{
					/* do linear search. */
					for ( l=cxadj[i]; l<cnedges; l++ )
					{
						if ( cadjncy[l] == k )
						{
							cadjwgt[l] += adjwgt[jj];
							break;
						}
					}	
					if ( l == cnedges )
					{
						cadjwgt[cnedges] = adjwgt[jj];
						cadjncy[cnedges] = k;
						htable[kk] = cnedges++;
					}
				}
			}
		}
		cxadj[i+1] = cnedges;

		if ( cxadj[i+1] - cxadj[i] < 1 )
		{
			printf("No neighbors for i:%d\n", i);
		}

		/* reset hashtable. Also compute cadjwgtsum. */
		cadjwgtsum[i] = 0;
		for ( j=cxadj[i]; j<cxadj[i+1]; j++ )
		{
			htable[cadjncy[j]&mask] = -1;
			cadjwgtsum[i] += cadjwgt[j];
		}
	}


	  /* Venu: my addition. Do not remove the contracted
	   * adjacency weight! Self loops are required for mlmcl to
	   * work. */
	  /* Remove the contracted adjacency weight */
	/*      jj = htable[cnvtxs&mask];
	  if (jj >= 0 && cadjncy[jj] != cnvtxs) {
		for (jj=0; jj<nedges; jj++) {
		  if (cadjncy[jj] == cnvtxs) 
			break;
		}
	  }
	  if (jj >= 0 && cadjncy[jj] == cnvtxs) {*/ /* This 2nd check is needed for non-adjacent matchings */
	/*        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
		cadjncy[jj] = cadjncy[--nedges];
		cadjwgt[jj] = cadjwgt[nedges];
	  } */

/*	ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));
*/

	cgraph->nedges = cnedges;

	ReAdjustMemory(graph, cgraph, dovsize);

	IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

	idxwspacefree(ctrl, mask+1);

}

/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);

  mask = HTLENGTH;
  if (cnvtxs < 8*mask || graph->nedges/graph->nvtxs > 15) { 
    CreateCoarseGraphNoMask(ctrl, graph, cnvtxs, match, perm);
    return;
  }

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
 /* cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
 */
  /* Venu: my addition */
  cgraph = my_SetUpCoarseGraph(graph, cnvtxs);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
/*    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      scopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);*/
    cvwgt[cnvtxs] = vwgt[v];

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m] += adjwgt[j];
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj] += adjwgt[j];
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = adjwgt[j];
        }
      }
    }

    if (v != u) { 
/*      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        saxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);
		*/
	  cvwgt[cnvtxs] += vwgt[u];

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m] += adjwgt[j];
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj] += adjwgt[j];
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = adjwgt[j];
          }
        }
      }

	  /* Venu: my addition. Do not remove the contracted
	   * adjacency weight! Self loops are required for mlmcl to
	   * work. */
      /* Remove the contracted adjacency weight */
/*      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) {*/ /* This 2nd check is needed for non-adjacent matchings */
/*        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      } */
    } 

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraphNoMask(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, k, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, dovsize;
  idxtype *xadj, *vwgt, *vsize, *adjncy, *adjwgt, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cvsize, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;

  dovsize = (ctrl->optype == OP_KVMETIS ? 1 : 0);
  dovsize=0; // Venu: my addtion here.

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  vsize = graph->vsize;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;


  /* Initialize the coarser graph */
 /* cgraph = SetUpCoarseGraph(graph, cnvtxs, dovsize);
 */
  /* Venu: my addition */
  cgraph = my_SetUpCoarseGraph(graph, cnvtxs);

  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cvsize = cgraph->vsize;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  htable = idxset(cnvtxs, -1, idxwspacemalloc(ctrl, cnvtxs));

  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    if (ncon == 1)
      cvwgt[cnvtxs] = vwgt[v];
    else
      scopy(ncon, nvwgt+v*ncon, cnvwgt+cnvtxs*ncon);
    cvwgt[cnvtxs] = vwgt[v];

    if (dovsize)
      cvsize[cnvtxs] = vsize[v];

    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      if ((m = htable[k]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = adjwgt[j];
        htable[k] = nedges++;
      }
      else {
        cadjwgt[m] += adjwgt[j];
      }
    }

    if (v != u) { 
      if (ncon == 1)
        cvwgt[cnvtxs] += vwgt[u];
      else
        saxpy(ncon, 1.0, nvwgt+u*ncon, 1, cnvwgt+cnvtxs*ncon, 1);

      if (dovsize)
        cvsize[cnvtxs] += vsize[u];

      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        if ((m = htable[k]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = adjwgt[j];
          htable[k] = nedges++;
        }
        else {
          cadjwgt[m] += adjwgt[j];
        }
      }

	  /* Venu: my addition. Do not remove the contracted
	   * adjacency weight! Self loops are required for mlmcl to
	   * work. */
	   /* But for some reason, removing the below code seems to
		* throw something into an infinite loop */
      /* Remove the contracted adjacency weight */
/*      if ((j = htable[cnvtxs]) != -1) {
        ASSERT(cadjncy[j] == cnvtxs);
        cadjwgtsum[cnvtxs] -= cadjwgt[j];
        cadjncy[j] = cadjncy[--nedges];
        cadjwgt[j] = cadjwgt[nedges];
        htable[cnvtxs] = -1;
      } */
    }

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d\n", cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt)));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]] = -1;  /* Zero out the htable */

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, dovsize);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, cnvtxs);
}


/*************************************************************************
* This function creates the coarser graph
**************************************************************************/
void CreateCoarseGraph_NVW(CtrlType *ctrl, GraphType *graph, int cnvtxs, idxtype *match, idxtype *perm)
{
  int i, j, jj, k, kk, l, m, istart, iend, nvtxs, nedges, ncon, cnedges, v, u, mask;
  idxtype *xadj, *adjncy, *adjwgtsum, *auxadj;
  idxtype *cmap, *htable;
  idxtype *cxadj, *cvwgt, *cadjncy, *cadjwgt, *cadjwgtsum;
  float *nvwgt, *cnvwgt;
  GraphType *cgraph;


  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->ContractTmr));

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  nvwgt = graph->nvwgt;
  adjncy = graph->adjncy;
  adjwgtsum = graph->adjwgtsum;
  cmap = graph->cmap;

  /* Initialize the coarser graph */
  cgraph = SetUpCoarseGraph(graph, cnvtxs, 0);
  cxadj = cgraph->xadj;
  cvwgt = cgraph->vwgt;
  cnvwgt = cgraph->nvwgt;
  cadjwgtsum = cgraph->adjwgtsum;
  cadjncy = cgraph->adjncy;
  cadjwgt = cgraph->adjwgt;


  iend = xadj[nvtxs];
  auxadj = ctrl->wspace.auxcore; 
  memcpy(auxadj, adjncy, iend*sizeof(idxtype)); 
  for (i=0; i<iend; i++)
    auxadj[i] = cmap[auxadj[i]];

  mask = HTLENGTH;
  htable = idxset(mask+1, -1, idxwspacemalloc(ctrl, mask+1)); 

  cxadj[0] = cnvtxs = cnedges = 0;
  for (i=0; i<nvtxs; i++) {
    v = perm[i];
    if (cmap[v] != cnvtxs) 
      continue;

    u = match[v];
    cvwgt[cnvtxs] = 1;
    cadjwgtsum[cnvtxs] = adjwgtsum[v];
    nedges = 0;

    istart = xadj[v];
    iend = xadj[v+1];
    for (j=istart; j<iend; j++) {
      k = auxadj[j];
      kk = k&mask;
      if ((m = htable[kk]) == -1) {
        cadjncy[nedges] = k;
        cadjwgt[nedges] = 1;
        htable[kk] = nedges++;
      }
      else if (cadjncy[m] == k) {
        cadjwgt[m]++;
      }
      else {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == k) {
            cadjwgt[jj]++;
            break;
          }
        }
        if (jj == nedges) {
          cadjncy[nedges] = k;
          cadjwgt[nedges++] = 1;
        }
      }
    }

    if (v != u) { 
      cvwgt[cnvtxs]++;
      cadjwgtsum[cnvtxs] += adjwgtsum[u];

      istart = xadj[u];
      iend = xadj[u+1];
      for (j=istart; j<iend; j++) {
        k = auxadj[j];
        kk = k&mask;
        if ((m = htable[kk]) == -1) {
          cadjncy[nedges] = k;
          cadjwgt[nedges] = 1;
          htable[kk] = nedges++;
        }
        else if (cadjncy[m] == k) {
          cadjwgt[m]++;
        }
        else {
          for (jj=0; jj<nedges; jj++) {
            if (cadjncy[jj] == k) {
              cadjwgt[jj]++;
              break;
            }
          }
          if (jj == nedges) {
            cadjncy[nedges] = k;
            cadjwgt[nedges++] = 1;
          }
        }
      }

	  /* Venu: my addition. Do not remove the contracted
	   * adjacency weight! Self loops are required for mlmcl to
	   * work. */
      /* Remove the contracted adjacency weight */
/*      jj = htable[cnvtxs&mask];
      if (jj >= 0 && cadjncy[jj] != cnvtxs) {
        for (jj=0; jj<nedges; jj++) {
          if (cadjncy[jj] == cnvtxs) 
            break;
        }
      }
      if (jj >= 0 && cadjncy[jj] == cnvtxs) {*/ /* This 2nd check is needed for non-adjacent matchings */
/*        cadjwgtsum[cnvtxs] -= cadjwgt[jj];
        cadjncy[jj] = cadjncy[--nedges];
        cadjwgt[jj] = cadjwgt[nedges];
      } */
    } 

    ASSERTP(cadjwgtsum[cnvtxs] == idxsum(nedges, cadjwgt), ("%d %d %d %d %d\n", cnvtxs, cadjwgtsum[cnvtxs], idxsum(nedges, cadjwgt), adjwgtsum[u], adjwgtsum[v]));

    for (j=0; j<nedges; j++)
      htable[cadjncy[j]&mask] = -1;  /* Zero out the htable */
    htable[cnvtxs&mask] = -1;

    cnedges += nedges;
    cxadj[++cnvtxs] = cnedges;
    cadjncy += nedges;
    cadjwgt += nedges;
  }

  cgraph->nedges = cnedges;

  ReAdjustMemory(graph, cgraph, 0);

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->ContractTmr));

  idxwspacefree(ctrl, mask+1);

}

/* Venu: my version of setupcoarsegraph*/
GraphType* my_SetUpCoarseGraph(GraphType *graph, int cnvtxs)
{
  GraphType *cgraph;

  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->ncon = graph->ncon;

  cgraph->finer = graph;
  graph->coarser = cgraph;
  /* make note of the fact that this is not the original graph */
  cgraph->isOrgGraph = 0;

  long size =sizeof(idxtype)*(3*cnvtxs+1 + 2*graph->nedges); 
  cgraph->gdata = (idxtype*)GKmalloc(size, "my_SetUpCoarseGraph: gdata");
//  cgraph->gdata = idxmalloc(3*cnvtxs+1 + 2*graph->nedges, "my_SetUpCoarseGraph: gdata");
  cgraph->xadj 	= cgraph->gdata;
  cgraph->adjwgtsum	= cgraph->gdata + cnvtxs+1;
  cgraph->cmap 	= cgraph->gdata + 2*cnvtxs+1;
//  cgraph->rmap2 = cgraph->gdata + 5*cnvtxs+1;
  cgraph->adjncy = cgraph->gdata + 3*cnvtxs+1;
  cgraph->adjwgt = cgraph->gdata + 3*cnvtxs+1 + graph->nedges;

  // we are allocating vwgt and rmap1 separately, as these two
  // are going to be needed for longer than gdata (i.e. we want
  // to free gdata, but don't want to free these two.
  cgraph->vwgt 	= idxmalloc(cnvtxs, "my_SetUpCoarseGraph: vwgt");
  cgraph->rmap1 = idxmalloc(cnvtxs, "my_SetUpCoarseGraph: rmap1");
 /* fprintf(stderr,"cgraph->rmap1:%d\n",(cgraph->rmap1==NULL?0:1));
  fprintf(stderr,"cgraph->rmap2:%d\n",(cgraph->rmap2==NULL?0:1));*/
  return cgraph;
}


/*************************************************************************
* Setup the various arrays for the coarse graph
**************************************************************************/
GraphType *SetUpCoarseGraph(GraphType *graph, int cnvtxs, int dovsize)
{
  GraphType *cgraph;

  cgraph = CreateGraph();
  cgraph->nvtxs = cnvtxs;
  cgraph->ncon = graph->ncon;

  cgraph->finer = graph;
  graph->coarser = cgraph;


  /* Allocate memory for the coarser graph */
  if (graph->ncon == 1) {
    if (dovsize) {
      cgraph->gdata = idxmalloc(5*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vwgt 		= cgraph->gdata + cnvtxs+1;
      cgraph->vsize 		= cgraph->gdata + 2*cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 3*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 5*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 5*cnvtxs+1 + graph->nedges;
    }
    else {
      cgraph->gdata = idxmalloc(4*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vwgt 		= cgraph->gdata + cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges;
    }
  }
  else {
    if (dovsize) {
      cgraph->gdata = idxmalloc(4*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->vsize 		= cgraph->gdata + cnvtxs+1;
      cgraph->adjwgtsum 	= cgraph->gdata + 2*cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 4*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 4*cnvtxs+1 + graph->nedges;
    }
    else {
      cgraph->gdata = idxmalloc(3*cnvtxs+1 + 2*graph->nedges, "SetUpCoarseGraph: gdata");
      cgraph->xadj 		= cgraph->gdata;
      cgraph->adjwgtsum 	= cgraph->gdata + cnvtxs+1;
      cgraph->cmap 		= cgraph->gdata + 2*cnvtxs+1;
      cgraph->adjncy 		= cgraph->gdata + 3*cnvtxs+1;
      cgraph->adjwgt 		= cgraph->gdata + 3*cnvtxs+1 + graph->nedges;
    }

    cgraph->nvwgt 	= fmalloc(graph->ncon*cnvtxs, "SetUpCoarseGraph: nvwgt");
  }

  return cgraph;
}


/*************************************************************************
* This function re-adjusts the amount of memory that was allocated if
* it will lead to significant savings
**************************************************************************/
void ReAdjustMemory(GraphType *graph, GraphType *cgraph, int dovsize) 
{

  if (cgraph->nedges > 100000 && graph->nedges < 0.7*graph->nedges) {
    idxcopy(cgraph->nedges, cgraph->adjwgt, cgraph->adjncy+cgraph->nedges);

    if (graph->ncon == 1) {
      if (dovsize) {
        cgraph->gdata = (idxtype*) realloc(cgraph->gdata, (5*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->vwgt 		= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->vsize 		= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 5*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 5*cgraph->nvtxs+1 + cgraph->nedges;
      }
      else {
        cgraph->gdata = (idxtype*) realloc(cgraph->gdata, (4*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 	= cgraph->gdata;
        cgraph->vwgt 	= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->cmap 	= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjncy 	= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjwgt 	= cgraph->gdata + 4*cgraph->nvtxs+1 + cgraph->nedges;
      }
    }
    else {
      if (dovsize) {
        cgraph->gdata = (idxtype*) realloc(cgraph->gdata, (4*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->vsize		= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->adjwgtsum	= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 4*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 4*cgraph->nvtxs+1 + cgraph->nedges;
      }
      else {
        cgraph->gdata = (idxtype*) realloc(cgraph->gdata, (3*cgraph->nvtxs+1 + 2*cgraph->nedges)*sizeof(idxtype));

        /* Do this, in case everything was copied into new space */
        cgraph->xadj 		= cgraph->gdata;
        cgraph->adjwgtsum	= cgraph->gdata + cgraph->nvtxs+1;
        cgraph->cmap 		= cgraph->gdata + 2*cgraph->nvtxs+1;
        cgraph->adjncy 		= cgraph->gdata + 3*cgraph->nvtxs+1;
        cgraph->adjwgt 		= cgraph->gdata + 3*cgraph->nvtxs+1 + cgraph->nedges;
      }
    }
  }

}

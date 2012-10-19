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
 * memory.c
 *
 * This file contains routines that deal with memory allocation
 *
 * Started 2/24/96
 * George
 *
 * $Id: memory.c,v 1.3 2009/12/07 06:01:33 venu Exp $
 *
 */

#include <metis.h>

void my_AllocateWorkSpace(CtrlType *ctrl, GraphType *graph )
{
  ctrl->wspace.pmat = NULL;

	long size = graph->nedges;
	size *= sizeof(EDegreeType);
//    ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(graph->nedges*sizeof(EDegreeType), "AllocateWorkSpace: edegrees");
    ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(size, "AllocateWorkSpace: edegrees");
    //ctrl->wspace.edegrees = (EDegreeType *)malloc(size);
	if ( ctrl->wspace.edegrees == NULL )
		errexit("Memory allocation failed for ctrl->wspace.edegrees. Requested size: %ld bytes", size);

	ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;
	ctrl->wspace.pmat = NULL;

	/* 4*nvtxs vectors for matching (see my_match_shemn function)
	 * and 5 nvtxs vectors for the intermediate rmaps */
//	 printf("nvtxs:%d\n",graph->nvtxs);
  	ctrl->wspace.maxcore = HTLENGTH+5*graph->nvtxs+10;
  	ctrl->wspace.core = idxmalloc(ctrl->wspace.maxcore, "AllocateWorkSpace: maxcore");
  ctrl->wspace.ccore = 0;
}

/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void AllocateWorkSpace(CtrlType *ctrl, GraphType *graph, int nparts)
{
  ctrl->wspace.pmat = NULL;

  printf("In AllocateWorkSpace\n");
  if (ctrl->optype == OP_KMETIS) {
    ctrl->wspace.edegrees = (EDegreeType *)GKmalloc(graph->nedges*sizeof(EDegreeType), "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases
          Coarsening
                    Matching: 4*nvtxs vectors
                 Contraction: 2*nvtxs vectors (from the above 4), 1*nparts, 1*Nedges
            Total = MAX(4*nvtxs, 2*nvtxs+nparts+nedges)

          Refinement
                Random Refinement/Balance: 5*nparts + 1*nvtxs + 2*nedges
                Greedy Refinement/Balance: 5*nparts + 2*nvtxs + 2*nedges + 1*PQueue(==Nvtxs)
            Total = 5*nparts + 3*nvtxs + 2*nedges

         Total = 5*nparts + 3*nvtxs + 2*nedges 
    */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           5*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else if (ctrl->optype == OP_KVMETIS) {
    ctrl->wspace.edegrees = NULL;
    ctrl->wspace.vedegrees = (VEDegreeType *)GKmalloc(graph->nedges*sizeof(VEDegreeType), "AllocateWorkSpace: vedegrees");
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.vedegrees;

    ctrl->wspace.pmat = idxmalloc(nparts*nparts, "AllocateWorkSpace: pmat");

    /* Memory requirements for different phases are identical to KMETIS */
    ctrl->wspace.maxcore = 3*(graph->nvtxs+1) +                 /* Match/Refinement vectors */
                           3*(nparts+1) +                       /* Partition weights etc */
                           graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* Greedy k-way balance/refine */
                           20  /* padding for 64 bit machines */
                           ;
  }
  else {
    ctrl->wspace.edegrees = (EDegreeType *)idxmalloc(graph->nedges, "AllocateWorkSpace: edegrees");
    ctrl->wspace.vedegrees = NULL;
    ctrl->wspace.auxcore = (idxtype *)ctrl->wspace.edegrees;

    ctrl->wspace.maxcore = 5*(graph->nvtxs+1) +                 /* Refinement vectors */
                           4*(nparts+1) +                       /* Partition weights etc */
                           2*graph->ncon*graph->nvtxs*(sizeof(ListNodeType)/sizeof(idxtype)) + /* 2-way refinement */
                           2*graph->ncon*(NEG_GAINSPAN+PLUS_GAINSPAN+1)*(sizeof(ListNodeType *)/sizeof(idxtype)) + /* 2-way refinement */
                           20  /* padding for 64 bit machines */
                           ;
  }

  ctrl->wspace.maxcore += HTLENGTH;
  ctrl->wspace.core = idxmalloc(ctrl->wspace.maxcore, "AllocateWorkSpace: maxcore");
  ctrl->wspace.ccore = 0;
}


/*************************************************************************
* This function allocates memory for the workspace
**************************************************************************/
void FreeWorkSpace(CtrlType *ctrl, GraphType *graph)
{
  GKfree( (void **)&ctrl->wspace.edegrees, (void **)&ctrl->wspace.vedegrees, (void **)&ctrl->wspace.core, (void **)&ctrl->wspace.pmat, LTERM);
}

/*************************************************************************
* This function returns how may words are left in the workspace
**************************************************************************/
int WspaceAvail(CtrlType *ctrl)
{
  return ctrl->wspace.maxcore - ctrl->wspace.ccore;
}

/*************************************************************************
* This function allocate space from the core 
**************************************************************************/
float *fwspacemalloc(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return (float *) (ctrl->wspace.core + ctrl->wspace.ccore - n);
}

/*************************************************************************
* This function frees space from the core 
**************************************************************************/
void fwspacefree(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}



/*************************************************************************
* This function creates a CoarseGraphType data structure and initializes
* the various fields
**************************************************************************/
GraphType *CreateGraph(void)
{
  GraphType *graph;

  graph = (GraphType *)GKmalloc(sizeof(GraphType), "CreateCoarseGraph: graph");

  InitGraph(graph);

  return graph;
}




void InitGraph(GraphType *graph) 
{
  graph->gdata = graph->rdata = NULL;

  graph->nvtxs = graph->nedges = -1;
  graph->mincut = graph->minvol = -1;

  graph->xadj = graph->vwgt = graph->adjncy = graph->adjwgt = NULL;
  graph->adjwgtsum = NULL;
  graph->label = NULL;
  graph->cmap = NULL;

  graph->where = graph->pwgts = NULL;
  graph->id = graph->ed = NULL;
  graph->bndptr = graph->bndind = NULL;
  graph->rinfo = NULL;
  graph->vrinfo = NULL;
  graph->nrinfo = NULL;

  graph->ncon = -1;
  graph->nvwgt = NULL;
  graph->npwgts = NULL;

  graph->vsize = NULL;

  graph->coarser = graph->finer = NULL;

  /* Venu: my addition */
  graph->rmap1=graph->rmap2=NULL;

}

void FreeGraph(GraphType *graph) 
{

  GKfree((void **)&graph->gdata, (void **)&graph->nvwgt, (void **)&graph->rdata, (void **)&graph->npwgts, LTERM);
  free(graph);
}

/*************************************************************************
* This function allocate space from the core 
**************************************************************************/
idxtype *idxwspacemalloc(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore += n;
  ASSERT(ctrl->wspace.ccore <= ctrl->wspace.maxcore);
  return ctrl->wspace.core + ctrl->wspace.ccore - n;
}

/*************************************************************************
* This function frees space from the core 
**************************************************************************/
void idxwspacefree(CtrlType *ctrl, int n)
{
  n += n%2; /* This is a fix for 64 bit machines that require 8-byte pointer allignment */

  ctrl->wspace.ccore -= n;
  ASSERT(ctrl->wspace.ccore >= 0);
}



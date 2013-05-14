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
 * match.c
 *
 * This file contains the code that computes matchings and creates the next
 * level coarse graph.
 *
 * Started 7/23/97
 * George
 *
 * $Id: match.c,v 1.2 2010/01/06 18:15:41 venu Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_RM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED && vwgt[i]+vwgt[adjncy[j]] <= ctrl->maxvwgt) {
          maxidx = adjncy[j];
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}

void my_Match_RM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt,*adjwgtsum;
  idxtype *match, *cmap, *perm, *rmap1, *rmap2;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum=graph->adjwgtsum;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  rmap1=idxwspacemalloc(ctrl,nvtxs);
  //rmap2=idxwspacemalloc(ctrl,nvtxs);
  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
  //      if (match[adjncy[j]] == UNMATCHED && vwgt[i]+vwgt[adjncy[j]] <= ctrl->maxvwgt) {
  		// Venu: my addition here
  		if (match[adjncy[j]] == UNMATCHED && adjncy[j] != i ) {
          maxidx = adjncy[j];
          break;
        }
      }

	  /* Venu: Added */
	  if ( maxidx == i )
	  {
	  	rmap1[cnvtxs]=i;
		//rmap2[cnvtxs]=-1;
	  }
	  else
	  {
		  if ( adjwgtsum[maxidx] > adjwgtsum[i] )
		  {
		  	rmap1[cnvtxs]=maxidx;
		  	//rmap2[cnvtxs]=i;
		  }
		  else
		  {
		  	rmap1[cnvtxs]=i;
		//	rmap2[cnvtxs]=maxidx;
	      }
	  }
	  /* End of Addition */

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  my_CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
//  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_RM_NVW(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *adjncy;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;

      /* Find a random matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED) {
          maxidx = adjncy[j];
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph_NVW(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}



/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_HEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, maxwgt;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *perm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  RandomPermute(nvtxs, perm, 1);

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
        if (match[k] == UNMATCHED && maxwgt < adjwgt[j] && vwgt[i]+vwgt[k] <= ctrl->maxvwgt) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}

void my_Match_PowerLaw_FC(CtrlType *ctrl, GraphType *graph)
{
	int i, ii, j, k, nvtxs, cnvtxs, maxidx, avgdegree,
		  numUnmatched,numMatched;
	idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum;
	idxtype *match, *cmap, *rmap1, *rmap2, *degrees, *perm, *tperm;

	nvtxs = graph->nvtxs;
	xadj = graph->xadj;
	vwgt = graph->vwgt; 
	adjncy = graph->adjncy;
	adjwgt = graph->adjwgt;
	adjwgtsum = graph->adjwgtsum;

	cmap = graph->cmap;
	rmap1 = idxset(nvtxs, -1, idxwspacemalloc(ctrl,nvtxs));
	match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));
	cmap = idxset(nvtxs, UNMATCHED, cmap);

	perm = idxwspacemalloc(ctrl, nvtxs);
	tperm = idxwspacemalloc(ctrl, nvtxs);
	degrees = idxwspacemalloc(ctrl, nvtxs);

	RandomPermute(nvtxs, tperm, 1); // Venu: my addition here.
	//    printf("Permuting according to degrees..\n");
	//	permuteDegreeOrder(nvtxs, tperm, graph->xadj);
	//   printf("Permuted according to degrees..\n");
	avgdegree = (int) 0.7*(xadj[nvtxs]/nvtxs);
	for (i=0; i<nvtxs; i++) 
	degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
	BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

	cnvtxs = 0;
	int reductionInSize = 0; 
	// This variable is for keeping track of how many
	// nodes have been eliminated during coarsening. We don't
	// want more than half the total number of nodes to be
	// eliminated in one coarsening step.

	idxtype* cvwgt = idxmalloc(nvtxs,
	"my_Match_PowerLaw_FC:cvwgt");
	// This one is for keeping track of how many nodes have been
	// merged together to form this current node. 

	idxtype* cpointers = idxmalloc(nvtxs,
	"my_Match_PowerLaw_FC:cpointers");
	// The constituent nodes of a parent node are stored in a
	// linked list using the match array, and this array points
	// to the first element of each list. i.e. cpointers[i] gives
	// the first element of the i^th parent node, and
	// match[cpointers[i]] gives the next element, and so on,
	// until -1 is reached.
	cpointers = idxset(nvtxs, UNMATCHED, cpointers);

	numUnmatched=0;
	numMatched=0;

	for ( ii=0 ; ii<nvtxs; ii++) 
	{
		i = perm[ii];

		if ( reductionInSize > nvtxs / 2  && cmap[i] == UNMATCHED)
		{
			cmap[i] = cnvtxs;
			cpointers[cnvtxs] = i;
			cvwgt[cnvtxs] = (vwgt == NULL) ? 1 : vwgt[i];
			rmap1[cnvtxs] = i;
			cnvtxs++;
			continue;
		}
	
		if (cmap[i] == UNMATCHED)   /* Unmatched */
		{
			maxidx = -1;
			float maxwgt = 0;
			//rtemp1 = 1.0/vwgt[i];
			float rtemp1 = 1.0/adjwgtsum[i];
			/* Find a heavy-edge matching, subject to maxvwgt constraints */
			for (j=xadj[i]; j<xadj[i+1]; j++) 
			{
				k = adjncy[j];
				//rtemp2 = adjwgt[j] *(rtemp1 + 1.0/vwgt[k]);
				float rtemp2;
				// rtemp2 = adjwgt[j] *(rtemp1 + 1.0/adjwgtsum[k]);

				// The following lines are a short-cut to sorting
				// first by edge wt followed by 1/vwgt. The 2nd
				// term should always be lesser than 1000, so it
				// should never make up for any difference in the
				// edge wts.
				if ( cmap[k] == UNMATCHED )
				{
					rtemp2 = 1000 * adjwgt[j] 
								+ 1.0 / (vwgt==NULL? 1 : vwgt[k]);
				}
				else
				{
					rtemp2 = 1000 * adjwgt[j] 
								+ 1.0 / cvwgt[cmap[k]];
				}

				if ( k != i && maxwgt < rtemp2 )// &&
					//(vwgt==NULL || vwgt[i]+vwgt[k] <= ctrl->maxvwgt) ) 
				{
					int valid = 1;
		/*			if ( cmap[k] == 517 )
					{
						printf("Hi there!\n");
					} */
					if ( cmap[k] != UNMATCHED )
					{
					if ( vwgt==NULL && 
							 cvwgt[cmap[k]]+1 > ctrl->maxvwgt )
							valid = 0;
						if ( vwgt != NULL && 
							 cvwgt[cmap[k]]+vwgt[i] >
							 ctrl->maxvwgt )
							 valid = 0;
					}
					else 
					{
						if ( vwgt != NULL && 
								vwgt[i]+vwgt[k] > ctrl->maxvwgt )
							valid = 0;
					}

					if ( valid )
					{
						maxwgt = rtemp2;
						maxidx = k;
					}
				}
			}

			if ( maxidx == -1)
			{
				cmap[i]=cnvtxs;
				rmap1[cnvtxs]=i;
				cpointers[cnvtxs] = i;
				if ( vwgt == NULL )
					cvwgt[cnvtxs] = 1;
				else
					cvwgt[cnvtxs] = vwgt[i];
				numUnmatched++;
				cnvtxs++;
			}
			else
			{
				reductionInSize++;
				if ( cmap[maxidx] == UNMATCHED )
				{
					// at this point both i and maxidx are unmatched.
					cmap[i] = cmap[maxidx] = cnvtxs;
					cpointers[cnvtxs] = i;
					match[i] = maxidx;
					cvwgt[cnvtxs] = (vwgt == NULL ) ? 2 :
										vwgt[i]	+ vwgt[maxidx];
					rmap1[cnvtxs] = i;
					cnvtxs++;
				}
				else
				{
					cmap[i] = cmap[maxidx];
					// put i at the head of the linked list
					// for cmap[maxidx].
					match[i] = cpointers[cmap[maxidx]];
					cpointers[cmap[maxidx]] = i;
					cvwgt[cmap[maxidx]] += (vwgt == NULL) ? 1 :
												vwgt[i] ;
					if ( cvwgt[cmap[maxidx]] > ctrl->maxvwgt )
					{
						printf("Oops! cvwgt[%d]:%d >",
								cmap[maxidx], cvwgt[cmap[maxidx]]);
						printf("maxvwgt:%d\n", ctrl->maxvwgt);
						fflush(stdout);
						abort();
					}
				}
			}
		}
		
/*		if ( ii == nvtxs-1)
		{
			printf("Done matching %d vertices.",ii);
			printf("NumMatched:%d, NumUnmatched:%d\n",numMatched,
			numUnmatched);
		} */
	}
	
	if ( reductionInSize + cnvtxs != nvtxs )
	{
		printf("reductionInSize (%d) + cnvtxs (%d) ", 
				reductionInSize, cnvtxs);
		printf("!= nvtxs (%d)\n", nvtxs);
		abort();
	}

	IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

	idxwspacefree(ctrl, nvtxs);  /* degrees */
	idxwspacefree(ctrl, nvtxs);  /* tperm */
	//printf("Done matching\n");

	CreateCoarseGraph_PowerLaw(ctrl, graph, cnvtxs, match,
			cpointers, cvwgt);
	/*
	fprintf(stderr,"graph->coarser:%d\n",(graph->coarser==NULL?0:1));
	fprintf(stderr,"graph->coarser->rmap1:%d\n",(graph->coarser->rmap1==NULL?0:1));
	fprintf(stderr,"graph->coarser->rmap2:%d\n",(graph->coarser->rmap2==NULL?0:1));
	*/
	memcpy(graph->coarser->rmap1, rmap1, cnvtxs*sizeof(idxtype));
	//  memcpy(graph->coarser->rmap2, rmap2, cnvtxs*sizeof(idxtype));

	free(cpointers);
	free(cvwgt);

	idxwspacefree(ctrl, nvtxs); /* rmap1 */
	//  idxwspacefree(ctrl, nvtxs); /* rmap2 */
	idxwspacefree(ctrl, nvtxs);
	idxwspacefree(ctrl, nvtxs);
}


/* Venu: my addition
*/
void my_Match_SHEMN(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, avgdegree,
		  numUnmatched,numMatched;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *match, *cmap, *rmap1, *rmap2, *degrees, *perm, *tperm;
  float rtemp1, rtemp2, maxwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

//  printf("Got to matching\n");

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt; 
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  cmap = graph->cmap;
  rmap1=idxwspacemalloc(ctrl,nvtxs);
  //rmap2=idxwspacemalloc(ctrl,nvtxs);
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);
//  degrees[nvtxs-1]=0;
// printf("xadj:%d\n",(xadj!=NULL));
//  printf("degrees:%d\n",(degrees!=NULL));
//  printf("nnz:%d\n",xadj[0]);

 // RandomPermute(nvtxs, tperm, 1);
  RandomPermute(nvtxs, tperm, 1); // Venu: my addition here.
//    printf("Permuting according to degrees..\n");
//	permuteDegreeOrder(nvtxs, tperm, graph->xadj);
 //   printf("Permuted according to degrees..\n");
  avgdegree = (int) 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) 
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  numUnmatched=0;
  numMatched=0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

	  /* Venu: Added */
	  if ( maxidx == i )
	  {
	  	rmap1[cnvtxs]=i;
		numUnmatched++;
		//rmap2[cnvtxs]=-1;
	  }
	  else
	  {
	  	  numMatched+=2;
		  if ( adjwgtsum[maxidx] > adjwgtsum[i] )
		  {
		  	rmap1[cnvtxs]=maxidx;
		 // 	rmap2[cnvtxs]=i;
		  }
		  else
		  {
		  	rmap1[cnvtxs]=i;
		//	rmap2[cnvtxs]=maxidx;
	      }
	  }
	  /* End of Addition */
	  
      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
/*	 if ( ii % 10000 == 0)
	 {
	 	printf("Done matching %d vertices.",ii);
		printf("NumMatched:%d, NumUnmatched:%d\n",numMatched,
		numUnmatched);
	  } */
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;
      //rtemp1 = 1.0/vwgt[i];
      rtemp1 = 1.0/adjwgtsum[i];
      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
	k = adjncy[j];
	//rtemp2 = adjwgt[j] *(rtemp1 + 1.0/vwgt[k]);
	rtemp2 = adjwgt[j] *(rtemp1 + 1.0/adjwgtsum[k]);
  //     if (match[k] == UNMATCHED && maxwgt < rtemp2 && vwgt[i]+vwgt[k] <= ctrl->maxvwgt) {
//       if (match[k] == UNMATCHED && maxwgt < rtemp2) {
	  // Venu: my addition here
        if (match[k] == UNMATCHED)
		{
		if ( k!=i && maxwgt < rtemp2) 
		{
		if (vwgt==NULL || vwgt[i]+vwgt[k] <= ctrl->maxvwgt)  
		{
          maxwgt = rtemp2;
          maxidx = adjncy[j];
        }
		}
		}
      }

	  /* Venu: Added */
	  if ( maxidx == i )
	  {
	  	rmap1[cnvtxs]=i;
		numUnmatched++;
		//rmap2[cnvtxs]=-1;
	  }
	  else
	  {
		  numMatched+=2;
		  if ( adjwgtsum[maxidx] > adjwgtsum[i] )
		  {
		  	rmap1[cnvtxs]=maxidx;
		  //	rmap2[cnvtxs]=i;
		  }
		  else
		  {
		  	rmap1[cnvtxs]=i;
		//	rmap2[cnvtxs]=maxidx;
	      }
	  }
	  /* End of Addition */
	  
      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
	 if ( ii == nvtxs-1)
	 {
//	 	printf("Done matching %d vertices.",ii);
//		printf("NumMatched:%d, NumUnmatched:%d\n",numMatched,
//		numUnmatched);
	  }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */
  //printf("Done matching\n");

  my_CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);
  /*
  fprintf(stderr,"graph->coarser:%d\n",(graph->coarser==NULL?0:1));
  fprintf(stderr,"graph->coarser->rmap1:%d\n",(graph->coarser->rmap1==NULL?0:1));
  fprintf(stderr,"graph->coarser->rmap2:%d\n",(graph->coarser->rmap2==NULL?0:1));
  */
  memcpy(graph->coarser->rmap1, rmap1, cnvtxs*sizeof(idxtype));
//  memcpy(graph->coarser->rmap2, rmap2, cnvtxs*sizeof(idxtype));

  idxwspacefree(ctrl, nvtxs); /* rmap1 */
//  idxwspacefree(ctrl, nvtxs); /* rmap2 */
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}

/* Venu: I am going to modify this function for my purposes */
void my_Match_HEMN(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt, *adjwgtsum;
  idxtype *match, *cmap, *rmap1, *rmap2, *perm;
  float rtemp1, rtemp2, maxwgt;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;
  adjwgtsum = graph->adjwgtsum;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));
 
  perm = idxwspacemalloc(ctrl, nvtxs);
  rmap1=idxwspacemalloc(ctrl,nvtxs);
//  rmap2=idxwspacemalloc(ctrl,nvtxs);
//  RandomPermute(nvtxs, perm, 1);
//  RandomPermute(nvtxs, perm, 1); // Venu: my addition here
	permuteDegreeOrder(nvtxs, perm, graph->xadj);
    printf("Permuted according to degrees..\n");

  cnvtxs = 0;
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;
      rtemp1 = 1.0/adjwgtsum[i];
      //rtemp1 = 1.0/vwgt[i];
      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        k = adjncy[j];
	rtemp2 = adjwgt[j] *(rtemp1 + 1.0/adjwgtsum[k]);
	//rtemp2 = adjwgt[j] *(rtemp1 + 1.0/vwgt[k]);
        if (match[k] == UNMATCHED && maxwgt < rtemp2 && vwgt[i]+vwgt[k] <= ctrl->maxvwgt) {
          maxwgt = rtemp2;
          maxidx = adjncy[j];
        }
      }

	  /* Venu: Added */
	  if ( maxidx == i )
	  {
	  	rmap1[cnvtxs]=i;
//		rmap2[cnvtxs]=-1;
	  }
	  else
	  {
		  if ( adjwgtsum[maxidx] > adjwgtsum[i] )
		  {
		  	rmap1[cnvtxs]=maxidx;
//		  	rmap2[cnvtxs]=i;
		  }
		  else
		  {
		  	rmap1[cnvtxs]=i;
//			rmap2[cnvtxs]=maxidx;
	      }
	  }
	  /* End of Addition */
	  
	  cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }
  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  my_CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);
  memcpy(graph->coarser->rmap1, rmap1, cnvtxs);
//  memcpy(graph->coarser->rmap2, rmap2, cnvtxs);

  idxwspacefree(ctrl, nvtxs);/* for rmap1*/
//  idxwspacefree(ctrl, nvtxs);/* for rmap2*/
  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


/*************************************************************************
* This function finds a matching using the HEM heuristic
**************************************************************************/
void Match_SHEM(CtrlType *ctrl, GraphType *graph)
{
  int i, ii, j, k, nvtxs, cnvtxs, maxidx, maxwgt, avgdegree;
  idxtype *xadj, *vwgt, *adjncy, *adjwgt;
  idxtype *match, *cmap, *degrees, *perm, *tperm;

  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->MatchTmr));

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  vwgt = graph->vwgt;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;

  cmap = graph->cmap;
  match = idxset(nvtxs, UNMATCHED, idxwspacemalloc(ctrl, nvtxs));

  perm = idxwspacemalloc(ctrl, nvtxs);
  tperm = idxwspacemalloc(ctrl, nvtxs);
  degrees = idxwspacemalloc(ctrl, nvtxs);

  RandomPermute(nvtxs, tperm, 1);
  avgdegree = (int) 0.7*(xadj[nvtxs]/nvtxs);
  for (i=0; i<nvtxs; i++) 
    degrees[i] = (xadj[i+1]-xadj[i] > avgdegree ? avgdegree : xadj[i+1]-xadj[i]);
  BucketSortKeysInc(nvtxs, avgdegree, degrees, tperm, perm);

  cnvtxs = 0;

  /* Take care any islands. Islands are matched with non-islands due to coarsening */
  for (ii=0; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      if (xadj[i] < xadj[i+1])
        break;

      maxidx = i;
      for (j=nvtxs-1; j>ii; j--) {
        k = perm[j];
        if (match[k] == UNMATCHED && xadj[k] < xadj[k+1]) {
          maxidx = k;
          break;
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  /* Continue with normal matching */
  for (; ii<nvtxs; ii++) {
    i = perm[ii];

    if (match[i] == UNMATCHED) {  /* Unmatched */
      maxidx = i;
      maxwgt = 0;

      /* Find a heavy-edge matching, subject to maxvwgt constraints */
      for (j=xadj[i]; j<xadj[i+1]; j++) {
        if (match[adjncy[j]] == UNMATCHED && maxwgt < adjwgt[j] && vwgt[i]+vwgt[adjncy[j]] <= ctrl->maxvwgt) {
          maxwgt = adjwgt[j];
          maxidx = adjncy[j];
        }
      }

      cmap[i] = cmap[maxidx] = cnvtxs++;
      match[i] = maxidx;
      match[maxidx] = i;
    }
  }

  IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->MatchTmr));

  idxwspacefree(ctrl, nvtxs);  /* degrees */
  idxwspacefree(ctrl, nvtxs);  /* tperm */

  CreateCoarseGraph(ctrl, graph, cnvtxs, match, perm);

  idxwspacefree(ctrl, nvtxs);
  idxwspacefree(ctrl, nvtxs);
}


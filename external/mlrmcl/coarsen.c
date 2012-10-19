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
 * coarsen.c
 *
 * This file contains the driving routines for the coarsening process 
 *
 * Started 7/23/97
 * George
 *
 * $Id: coarsen.c,v 1.3 2010/03/12 00:25:53 venu Exp $
 *
 */

#include <metis.h>


/*************************************************************************
* This function takes a graph and creates a sequence of coarser graphs
**************************************************************************/
GraphType *Coarsen2Way(CtrlType *ctrl, GraphType *graph)
{
  int clevel;
  GraphType *cgraph;

//  IFSET(ctrl->dbglvl, DBG_TIME, starttimer(ctrl->CoarsenTmr));
  cleartimer(ctrl->CoarsenTmr);
  starttimer(ctrl->CoarsenTmr);

  cgraph = graph;

  if ( ctrl->CoarsenTo >= graph->nvtxs )
  {
  	// don't do any coarsening.
	return cgraph;
  }

  /* The following is ahack to allow the multiple bisections to go through with correct
     coarsening */
  if (ctrl->CType > 20) {
    clevel = 1;
    ctrl->CType -= 20;
  }
  else
    clevel = 0;

  do {
    IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d [%d] [%d %d]\n",
          cgraph->nvtxs, cgraph->nedges, ctrl->CoarsenTo, ctrl->maxvwgt, 
          (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt) : cgraph->nvtxs)));

    if (cgraph->adjwgt) {
      switch (ctrl->CType) {
        case MATCH_RM:
          my_Match_RM(ctrl, cgraph);
          break;
        case MATCH_HEM:
          if (clevel < 1)
            Match_RM(ctrl, cgraph);
          else
            Match_HEM(ctrl, cgraph);
          break;
        case MATCH_HEMN:
          if (clevel < 1)
            Match_RM(ctrl, cgraph);
          else
            my_Match_HEMN(ctrl, cgraph);
          break;
        case MATCH_SHEMN:
          //if (clevel < 1)
          //  Match_RM(ctrl, cgraph); 
          //else
            my_Match_SHEMN(ctrl, cgraph);
          break;
		case MATCH_POWERLAW_FC:
			my_Match_PowerLaw_FC(ctrl, cgraph);
		  break;
        case MATCH_SHEM:
          if (clevel < 1)
            Match_RM(ctrl, cgraph);
          else
            Match_SHEM(ctrl, cgraph);
          break;
        case MATCH_SHEMKWAY:
          Match_SHEM(ctrl, cgraph);
          break;
/*		case MATCH_HASH:
			printf("clevel:%d\n",clevel);
			if ( clevel < 3 )	
		 		match_hash(ctrl, cgraph,4);
			else
				my_Match_SHEMN(ctrl, cgraph);
		  break; */
        default:
          errexit("Unknown CType: %d\n", ctrl->CType);
      }
    }
    else {
      Match_RM_NVW(ctrl, cgraph);
    }

    cgraph = cgraph->coarser;
    clevel++;
//	printf("cgraph->nvtxs::%d, ", cgraph->nvtxs);
//	printf("cgraph->nedges::%d\n", cgraph->nedges);
//	fflush(stdout);

	/* Venu: my addition */
/*  } while (cgraph->nvtxs > ctrl->CoarsenTo); 
*/
  } while ((cgraph->nvtxs > ctrl->CoarsenTo && cgraph->nvtxs <
  COARSEN_FRACTION2*cgraph->finer->nvtxs && cgraph->nedges >
  cgraph->nvtxs/2) );
  			//|| cgraph->nvtxs > 50000 ); 

  if ( cgraph->nvtxs > ctrl->CoarsenTo )
  {
  	if ( cgraph->nedges <= cgraph->nvtxs/2 )
	{
		printf("Stopped coarsening as nedges (%d) is ", cgraph->nedges );
		printf("fewer than half of nvtxs (%d)\n", cgraph->nvtxs	);
	}
	else if ( cgraph->nvtxs > COARSEN_FRACTION2*cgraph->finer->nvtxs )
	{
		printf("Stopped coarsening as ratio of new nvtxs");
		printf(" (%d) to old nvtxs (%d) > %.2f\n", 
				cgraph->nvtxs, cgraph->finer->nvtxs,
				COARSEN_FRACTION2 );
	}
  }


  IFSET(ctrl->dbglvl, DBG_COARSEN, printf("%6d %7d [%d] [%d %d]\n",
        cgraph->nvtxs, cgraph->nedges, ctrl->CoarsenTo, ctrl->maxvwgt, 
        (cgraph->vwgt ? idxsum(cgraph->nvtxs, cgraph->vwgt) : cgraph->nvtxs)));

  //IFSET(ctrl->dbglvl, DBG_TIME, stoptimer(ctrl->CoarsenTmr));
  stoptimer(ctrl->CoarsenTmr);

 // if ( ctrl->CType == MATCH_HASH )
 // 	freePrimes();

  return cgraph;
}


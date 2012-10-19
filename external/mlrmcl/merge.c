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

void addToList(ListInt *l, int a)
{
	if ( l->length+1 > l->allocSize )
	{
		l->allocSize += l->increment;
		l->l = (idxtype*)realloc(l->l, sizeof(idxtype)*(l->allocSize));
	}
	l->l[l->length++] = a;
}

void addToList(ListWgt *l, int a)
{
	if ( l->length+1 > l->allocSize )
	{
		l->allocSize += l->increment;
		l->l = (wgttype*)realloc(l->l, sizeof(wgttype)*(l->allocSize));
	}
	l->l[l->length++] = a;
}

idxtype* convertMetisOutputToClusterLists(const idxtype* part, int
		nvtxs, int nparts, idxtype** part_xadj)
{
	*part_xadj = idxmalloc(nparts+1, "part_xadj");

	for ( int i = 0; i<nparts+1; i++ )
		(*part_xadj)[i] = 0;

	for ( int i = 0; i<nvtxs; i++ )
		(*part_xadj)[part[i]+1]++;

	(*part_xadj)[0]=0;
	for ( int i=1; i<nparts+1; i++ )
		(*part_xadj)[i] = (*part_xadj)[i-1] + (*part_xadj)[i];
	
	idxtype *part_rev = idxmalloc(nvtxs+1, "part_rev");
	idxtype *counters = idxmalloc(nparts, "counters");
	for ( int i=0; i<nparts; i++ )
		counters[i] = 0;

	for ( int i=0; i<nvtxs; i++ )
	{
		int pid = part[i];
		part_rev[counters[pid] + (*part_xadj)[pid]] = i;
		counters[pid]++;
	}

	free(counters);
	return part_rev;
}

void copyListWgt(ListWgt *dest, ListWgt *src)
{
	dest->l = src->l;
	dest->length = src->length;
	dest->allocSize = src->allocSize;
	dest->increment = src->increment;
}

void copyListInt(ListInt *dest, ListInt *src)
{
	dest->l = src->l;
	dest->length = src->length;
	dest->allocSize = src->allocSize;
	dest->increment = src->increment;
}

void allocListWgt(ListWgt *l, int initSize, int increment)
{
	l->l = (wgttype*) malloc(sizeof(wgttype)*initSize);
	l->length = 0;
	l->allocSize = initSize;
	l->increment = increment;
}

void allocListInt(ListInt *l, int initSize, int increment)
{
	l->l = idxmalloc(initSize, "allocListInt");
	l->length = 0;
	l->allocSize = initSize;
	l->increment = increment;
}

ListGraph* allocListGraph(int nvtxs, int wgtflag, int
initAdjSize, int increment)
{
	ListGraph* lg= (ListGraph*)malloc(sizeof(ListGraph));
	lg->nvtxs = nvtxs;
	lg->adjLists = (ListInt*)malloc(sizeof(ListInt)*nvtxs);
	if ( wgtflag )
	{
		lg->wgtLists = (ListWgt*)malloc(sizeof(ListWgt)*nvtxs);
	}

	for ( int i=0; i<nvtxs; i++ )
	{
		allocListInt(&(lg->adjLists[i]), initAdjSize, increment);

		if ( wgtflag )
		{
			allocListWgt(&(lg->wgtLists[i]), initAdjSize,
			increment);
		}
	}

	lg->vols = (wgttype*) malloc(sizeof(wgttype) * nvtxs);

	return lg;
}

ListGraph* createClusterGraph(const idxtype *part, int nparts,
				const GraphType* g)
{
	ListGraph* cg = allocListGraph(nparts, 1, 100, 100);
	
	idxtype *partLists, *part_xadj;
	partLists = convertMetisOutputToClusterLists(part, g->nvtxs,
				nparts, &part_xadj);

	wgttype* dense_wgts = (wgttype*)malloc(sizeof(wgttype)*nparts);
	int totalEdges = 0;
	for ( int pid = 0; pid < nparts; pid++)
	{
		cg->vols[pid] = 0;
		// initialize dense wgt vector
		for ( int i = 0; i < nparts; i++ )
			dense_wgts[i] = 0;

	   // determine adjacencies of pid, and the edge weights
		for ( int pid_j = part_xadj[pid]; pid_j < part_xadj[pid+1]; pid_j++ )
		{
			int member = partLists[pid_j];
			for ( int j = g->xadj[member]; j < g->xadj[member+1]; j++ )
			{
				int nbr_part = part[g->adjncy[j]];
				wgttype nbr_wgt = 1;
				if ( g->adjwgt != NULL )
					nbr_wgt = g->adjwgt[j];
				
				//if ( nbr_part == pid )
				//	nbr_wgt /= 2; 
					// since these edges will be repeating in our
					// calculations of vols and selfassoc.
			
				cg->vols[pid] += nbr_wgt; 
				// total volume of cluster
			 	
				// only store weights in the upper triangular part,
				// since the graph's symmetric anyway
				if ( nbr_part >= pid )
				{
					 dense_wgts[nbr_part] += nbr_wgt; 
				}
			}
		}

		// now store in sparse format
		for ( int i = 0; i < nparts; i++ )
		{
			if ( dense_wgts[i] > 0 || i == pid )
			{
				addToList(&(cg->adjLists[pid]), i);
				addToList(&(cg->wgtLists[pid]), dense_wgts[i]);
				totalEdges++;
			}
		}
		if ( (cg->adjLists[pid]).l[0] != pid )
		{
			printf("Yikes! First neighbor is not self: %d\n",
			pid);
			abort();
		}
	}
	free(dense_wgts);
	cg->nedges = totalEdges;

	// now normalize each edge weight by the cluster volumes
/*	for ( int pid = 0; pid < nparts; pid++ )
	{
		for ( int j = 0; j < cg->adjLists[pid]->length; j++ )
		{
			int nbr = cg->adjLists[pid]->l[j];
			wgttype maxVol = (vols[nbr] > vols[pid]) ? vols[nbr] :
									vols[pid];
			cg->wgtLists[pid]->l[j] /= maxVol;
		}
	} */

	free( partLists );

	return cg;
}

ListInt* mergeLists(const ListInt* al1, const ListWgt* wl1, 
		const ListInt* al2, const ListWgt* wl2, ListWgt** ret_wl,
		int univ_size)
{
	ListInt* ret_al = (ListInt*) malloc(sizeof(ListInt));
	*ret_wl = (ListWgt*) malloc(sizeof(ListWgt));

	allocListInt(ret_al, al1->length+al2->length, 100);
	allocListWgt(*ret_wl, al1->length+al2->length, 100);

	// we need to have hashtables
	idxtype* ht1, *ht2;
	ht1 = idxmalloc(univ_size, "ht1");
	ht2 = idxmalloc(univ_size, "ht2");
	for ( int i=0; i<univ_size; i++ )
		ht1[i] = ht2[i] = -1;

	for ( int i=0; i<al1->length; i++ )
		ht1[al1->l[i]] = i;

	for ( int i=0; i<al2->length; i++ )
		ht2[al2->l[i]] = i;

	// go through al1 first, making sure to add contributions
	// from al2
	for ( int i=0; i<al1->length; i++ )
	{
		wgttype wt = wl1->l[i];
		int t;
		if ( (t = ht2[al1->l[i]]) > -1 )
			wt += wl2->l[t];
		addToList(ret_al, al1->l[i]);
		addToList(*ret_wl, wt);
	}

	// go through al2 now, and *only* include elements that were
	// absent in al1, since the others have already been taken
	// into account.
	for ( int i=0; i<al2->length; i++ )
	{
		if ( ht1[al2->l[i]] > -1 )
			continue;
		addToList(ret_al, al2->l[i]);
		addToList(*ret_wl, wl2->l[i]);
	}
	
	free(ht1);
	free(ht2);

	return ret_al;
}

void mergeBestClusters(ListGraph *cg, idxtype* part_mapper, int
									nparts, int mergeHeuristic)
{
	wgttype maxScore = 0;
	idxtype pid1=-1, pid2=-1;

	// go through all edge weights and figure out the best.
	for ( int i=0; i<nparts; i++)
	{
		if ( part_mapper[i] != i )
			continue;
		wgttype selfassoc = (cg->wgtLists[i]).l[0];
		for ( int j=0; j<(cg->wgtLists[i]).length; j++ )
		{
			int nbr = (cg->adjLists[i]).l[j];
			if ( part_mapper[nbr] != nbr || nbr == i )
				continue;

			wgttype score;
			if (  mergeHeuristic == 2 )
			{
				 wgttype cross = (cg->wgtLists[i]).l[j];
				 wgttype nbrassoc = (cg->wgtLists[nbr]).l[0];
				 wgttype score = 
				 (selfassoc + nbrassoc + cross)/(cg->vols[i] + cg->vols[nbr])
				 - selfassoc/cg->vols[i] - nbrassoc/cg->vols[nbr];
			}
			else
			{
				 wgttype maxVol = (cg->vols[i] < cg->vols[nbr]) ?
										 cg->vols[nbr] : cg->vols[i] ;
				 score = (cg->wgtLists[i]).l[j] / maxVol; 
			}

			if ( pid1 == -1 || score > maxScore )
			{
				maxScore = score;
				pid1 = i;
				pid2 = nbr;
			} 
		}
	}

	if ( pid1 == -1 )
	{
		printf("Yikes! All clusters seem to be disconnected\n");
		return;
	}

	printf("Going to merge %d and %d clusters\n", (pid1+1),
				(pid2+1));

	// merge adj lists
	ListWgt* new_wl; 
	ListInt* new_al = mergeLists(&cg->adjLists[pid1],
	&cg->wgtLists[pid1], &cg->adjLists[pid2], &cg->wgtLists[pid2],
	&new_wl, nparts);

	// pid1 and pid2 will be mapped to smaller of pid1 and pid2;
	int minpid = (pid1 < pid2) ? pid1 : pid2;
	int maxpid = (pid1 == minpid) ? pid2 : pid1;
	part_mapper[pid1] = part_mapper[pid2] = minpid;
	cg->vols[minpid] = cg->vols[pid1] + cg->vols[pid2];
	// for all partitions mapped to maxpid, update them to minpid
	for ( int i=0; i<nparts; i++ )
	{
		if ( part_mapper[i] == maxpid )
			part_mapper[i] = minpid;
	}	

	// dealloc existing minpid lists
	free((cg->adjLists[minpid]).l);
//	free(cg->adjLists[minpid]);
//	cg->adjLists[minpid] = NULL;
	free((cg->wgtLists[minpid]).l);
//	free(cg->wgtLists[minpid]);
//	cg->wgtLists[minpid] = NULL;

	// assign new lists..
	copyListInt(&(cg->adjLists[minpid]), new_al);
	copyListWgt(&(cg->wgtLists[minpid]), new_wl);
//	cg->adjLists[minpid]->l = new_al->l;
//	cg->wgtLists[minpid]->l = new_wl->l;

	// the edges between clusters with id less than minpid will
	// not have been update in new_al; we'll have to update them
	// separately (since we're only storing upper triangular
	// matrix)
	for ( int i=0; i<minpid; i++ )
	{
		wgttype newwt = 0;
		idxtype minpidid = -1;
		for ( int j=0; j<(cg->adjLists[i]).length; j++ )
		{
			if ( (cg->adjLists[i]).l[j] == minpid ||
					(cg->adjLists[i]).l[j] == maxpid )
			{
				newwt += (cg->wgtLists[i]).l[j];
				if ( (cg->adjLists[i]).l[j] == minpid )
					minpidid = j;
				else
					break; 
			}
		}
		if ( minpidid > -1 )
			(cg->wgtLists[i]).l[minpidid] = newwt;
	}

}

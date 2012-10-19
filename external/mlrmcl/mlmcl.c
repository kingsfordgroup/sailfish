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


/* This file has the top-level routines for multi-level MCL
 * algorithm 
 * 
 */
 
#include <metis.h>

/**
 * NOTE!! THIS FUNCTION IS NO LONGER BEING USED, IT HAS BEEN
 * SUPERSEDED BY dprmcl in rmcl.c
 **/
/*
Matrix* mclForMultiLevel(Matrix* flows, Matrix* adj, float
gamma, int maxiter, int dump, int exact, int k, int level,wgttype
dpr_threshold, GraphType* graph)
{
	Matrix *M0,*M, *Mnew;
	int numiter,i;
	int diff;
	int numUnique=0;
	int numWithoutAttrs=0;
	int calcChange=0,getNcut=0,dprPhase=0,calcDiff=0;
	wgttype changeInNorm=-1,changeInNcut=-1,prevNcut=-1;
	timer expTimer,infTimer,pnTimer,totalTimer;
	idxtype* prevAttractors=idxmalloc(adj->nvtxs,"Attractors");
	idxtype* hashtable=idxmalloc(adj->nvtxs,"Hashtable");

	// clear hashtable.
	for(i=0;i<adj->nvtxs;i++)
		hashtable[i]=-1;
	
	M0=adj;
	M=flows;

	numiter=0;
	do
	{
		if ( dprPhase <= 0 && numiter > 9 && level <= 0 
			&& (changeInNorm != -1 && changeInNorm < 1.0e-03) )
		//	&& (changeInNcut != -1 && changeInNcut < 5.0e-03) )
		{
			dprPhase=1;
		}

		if ( numiter > 9 && level <= 0 )
		{
			calcChange=1; //getNcut=1;
		}

		cleartimer(totalTimer);

		starttimer(totalTimer);

		if ( dprPhase > 0 )
		{
			M0=getDprAdjMatrix(M,adj,hashtable,dpr_threshold);
			printf("nnz of actual adj matrix:%d,",adj->nnz);
			printf("nnz of dpr adj matrix:%d\n",M0->nnz);
		} 

		Mnew=allInOneStep(M,M0,hashtable,gamma,exact,k,calcChange,0);

		if ( calcChange > 0 )
		{
			changeBetweenMatrices(Mnew,M,&changeInNorm);///(M->nvtxs);
			changeInNorm=changeInNorm/M->nvtxs;
		}

		stoptimer(totalTimer);

		if ( M != M0 )
			freeMatrix(M);
		if ( M0 != NULL &&  M0 != adj )
		{
			freeMatrix(M0);
			M0=adj;
		}
		

//		if ( dump > 0 )
//			dumpMatrix(Mnew);
		M=Mnew;

		if ( calcDiff > 0 )
		{
//			numUnique=noOfUnique(Mnew->attractors,Mnew->nvtxs,hashtable);
			if ( numiter > 0 )
				diff=isConverged2(prevAttractors,M->attractors,M->nvtxs);
			else	
				diff=M->nvtxs;
			numWithoutAttrs=numWithoutAttractors(M);
		}

		if ( getNcut > 0 )
		{
			float ncut=0;
		//	idxtype* attrs = (idxtype*)
		//			malloc(sizeof(idxtype)*M->nvtxs);
			idxtype* attrs=idxmalloc(M->nvtxs,
							"mclForMultilevel:attrs for ncut");
			getAttractorsForAll(M);
			idxcopy(M->nvtxs,M->attractors,attrs);
			numUnique=mapPartition(attrs,M->nvtxs);
			ncut=ComputeNCut(graph,attrs,numUnique)/numUnique;
			printf("num. clusters:%d, avg. ncut:%.3f\n",
						numUnique,ncut);
			if ( prevNcut != -1 )
			{
				changeInNcut = (prevNcut - ncut)/(prevNcut);
			}
			free(attrs);
			prevNcut=ncut;
		}
	
		idxcopy(M->nvtxs,M->attractors,prevAttractors);

		if ( numiter%1 == 0 )
		{
			if ( calcDiff > 0 )
			{
				printf("diff:%d,numWithoutAttrs:%d\n",
					diff,numWithoutAttrs);
			}
			if ( calcChange > 0 )
				printf("avg. change:%.6f\n",changeInNorm);
		
			printf("Done with %d iterations\n", numiter);
			printf("Total time for iteration:%7.3f\n",
						gettimer(totalTimer)); 

			fflush(stdout);
		} 
		numiter++;

		if ( calcChange > 0 && changeInNorm < 1.0e-05 )
			break;

	}while(	numiter < maxiter );
//	}while( !isConverged(prevAttractors,M->attractors,M->nvtxs)
//				&& (numiter < maxiter) );
//	}while(diff>(0.001)*((double)M->nvtxs) && numiter<maxiter);
//	}while( (numWithoutAttrs>(0.0001*((double)M->nvtxs)) 
//			|| diff>(0.0001)*((double)M->nvtxs)) && numiter<maxiter);

	
//	numUnique=noOfUnique(Mnew->attractors,Mnew->nvtxs,hashtable);
//	printf("no. of attractors:%d\n", numUnique);
	
	free(hashtable);
	free(prevAttractors);
	
	// This is just to verify that M0 has not been freed. If M0
	// has been freed, this will throw an exception. 
//	printf("nnz:%d\n", M0->nnz); 

	return M;
}
*/

Matrix* propagateFlow(Matrix* cm, GraphType* cgraph, GraphType*
rgraph, int nnzRefinedGraph)
{
	int i,j,k,nnzcount=0;
	idxtype* cmap=rgraph->cmap,*rmap1=cgraph->rmap1,
				*rmap2=cgraph->rmap2;

	Matrix* rm=allocMatrix(rgraph->nvtxs,nnzRefinedGraph,0,0,0);

	rm->xadj[0]=0;
	for(i=0;i<rm->nvtxs;i++)
	{
		j=cmap[i];
		for(k=cm->xadj[j];k<cm->xadj[j+1];k++)
		{
			if ( rmap1[cm->adjncy[k]] >= rgraph->nvtxs )
			{
				printf("Oops! rmap1[%d] = %d,", cm->adjncy[k],
							rmap1[cm->adjncy[k]] );
				printf(" > %d\n", rgraph->nvtxs);
				abort();
			}
			if ( nnzcount >= nnzRefinedGraph )
			{
				printf("Incorrect calculation of");
				printf(" nnzRefinedGraph in propagateFlow\n");
				printf("nodeId:%d, rm->nvtxs:%d, ", i,
						rm->nvtxs);
				printf("nnzRefinedGraph:%d\n", nnzRefinedGraph);
				abort();
			}
			rm->adjncy[nnzcount]=rmap1[cm->adjncy[k]];
			rm->adjwgt[nnzcount]=cm->adjwgt[k];
			nnzcount++;
		}
		rm->xadj[i+1]=nnzcount;
	}
	rm->nnz=nnzcount;

//	printf("Propagated %d flow values\n",nnzcount);
//	printf("avg nnz per column:%d\n",(nnzcount/rm->nvtxs));

	return rm;
}

void mlmclWithGraph(GraphType *graph,idxtype* indices,Options opt)
{
	GraphType *cgraph,*t;
	int nvtxs=graph->nvtxs;
	CtrlType ctrl;
	int levels=0;
	Matrix *M=NULL;

	int ct = ctrl.CoarsenTo = opt.coarsenTo, runMcl=0;
	ctrl.CType=opt.matchType;
/*	if ( nvtxs > 50000 )
		ctrl.CType = MATCH_POWERLAW_FC; */
	ctrl.optype=OP_KMETIS;
//	ctrl.maxvwgt = (int) round( ((ct>100)?nvtxs/100: nvtxs/ct));
//	ctrl.maxvwgt = (ctrl.maxvwgt > 5000 ) ? 5000 : ctrl.maxvwgt ;
	ctrl.dbglvl=1;

	if ( ct < 5 )
		ctrl.maxvwgt = (int) ceil( (1.5*nvtxs)/ct );
	else if ( ct <= 100 )
		ctrl.maxvwgt = (int) ceil( (2.0*nvtxs)/ct );
	else if ( ct > 100 )
	{
		// we can allow some imbalance here.
		ctrl.maxvwgt = (int) ceil( 10.0 * nvtxs/ct );
	}

	my_AllocateWorkSpace(&ctrl,graph);

	cgraph = Coarsen2Way(&ctrl,graph);	

	FreeWorkSpace(&ctrl, graph);

//	printf( "Coarsen time:%.2f\n",gettimer(ctrl.CoarsenTmr) );
	for(levels=0,t=cgraph;t->finer!=NULL;t=t->finer,levels++);

	printf("Coarsened %d levels\n", levels);
	printf("Number of vertices in coarsest graph:%d\n", cgraph->nvtxs);
//	printf("In lib/mlmcl.c: penalty_power:%.1f\n", opt.penalty_power);
	fflush(stdout);

	while ( cgraph != NULL )
	{
		Matrix *M0, *Mnew;
		int iters;
//		dump_graph(cgraph);
		M0=setupCanonicalMatrix(cgraph->nvtxs, cgraph->nedges,
			cgraph->xadj, cgraph->adjncy, cgraph->adjwgt, opt.ncutify);

		if ( cgraph->coarser != NULL )
		{
			// if this is not the original graph, then we don't
			// need the adjacency information in the graph
			// anymore. M0 is what we want, and we have it, so we
			// can free cgraph->gdata.
			GKfree( (void**)&(cgraph->gdata), LTERM);
		}

		if ( cgraph->coarser == NULL )
		{
		// if this is the coarsest graph, flow matrix is same as
		// the transition matrix. 
			M=M0;
	//		dumpMatrix(M);
		}

		if (cgraph->finer == NULL)
			iters=opt.num_last_iter;
		else
			iters=opt.iter_per_level;

		if ( cgraph->coarser == NULL )
		{;
			//printRow(M,35);
			//printRow(M,36);
		}

		Mnew=dprmcl(M,M0,cgraph, opt, iters, levels);
		// M would have been freed in mclForMultiLevel, so no
		// need to free it here.
	/*	if ( M0 != NULL )
		{
			freeMatrix(M0);
		} */

		if ( cgraph->coarser == NULL )
		{
			//printRow(Mnew,35);
			//printRow(Mnew,36);
		}

		if ( runMcl > 0 && cgraph->finer == NULL )
		{
			/* run original mcl so that the matrix moves closer to
			 * convergence */
			M=dprmcl(Mnew,Mnew,cgraph, opt, iters/2,levels);
			/* Mnew would not have been freed in mclForMultiLevel
			 * above. */
			freeMatrix(Mnew);
		}
		else
		{
			M=Mnew;
		}

		if ( cgraph->finer != NULL )
		{
			int nnzRefinedGraph = 2*M->nnz;
			if ( ctrl.CType == MATCH_POWERLAW_FC )
			{
				int ii;
				nnzRefinedGraph = 0;
				for ( ii=0; ii<cgraph->finer->nvtxs; ii++ )
				{
					int tx=cgraph->finer->cmap[ii];
					nnzRefinedGraph +=
					M->xadj[tx+1]-M->xadj[tx];
				}
			}
			else if ( ctrl.CType == MATCH_HASH  && levels <= 3 )
			{
				int ii;
				nnzRefinedGraph=0;
				for ( ii=0; ii<cgraph->nvtxs; ii++ )
				{
					nnzRefinedGraph += cgraph->vwgt[ii] *
							(M->xadj[ii+1] - M->xadj[ii]);
				}
			}
//			dumpMatrix(M);
			Mnew=propagateFlow(M,cgraph,cgraph->finer,nnzRefinedGraph);
//			dumpMatrix(Mnew);
			if ( M != NULL )
			{
				freeMatrix(M);
			}
			M=Mnew;
		}
		cgraph=cgraph->finer;

		// These two didn't get freed earlier, when we freed
		// gdata.
		if ( cgraph != NULL )
		{
			GKfree( (void**)&(cgraph->coarser->vwgt),
					(void**)&(cgraph->coarser->rmap1), LTERM);
		}

		levels--;
		printf("Done level %d\n", levels+1);
		fflush(stdout);
	}

	getAttractorsForAll(M);
	idxcopy(nvtxs, M->attractors, indices);

/*	int nClusters = 0;
	wgttype minWgt = 0;
	idxtype *ret = getNodesToComponentMap(M, &nClusters, minWgt);
	idxcopy(nvtxs, ret, indices);
	free(ret);
*/
	freeMatrix(M);
}

Matrix* mis_projectFlow(Matrix* Mc, idxtype* rmap, int nvtxs)
{
	Matrix* Mr = allocMatrix(nvtxs, Mc->nnz, 0, 0, 0);
	idxtype* hashtable = idxmalloc(nvtxs,
						"mis_projectFlow:hashtable");
	int i,Mr_counter;
	
	for ( i=0; i<nvtxs; i++ )
	{
		hashtable[i]=-1;
	}

	for ( i=0; i<Mc->nvtxs; i++ )
	{
		hashtable[rmap[i]]=i;
	}

	Mr->xadj[0]=0;
	Mr_counter=0;
	for ( i=0; i<nvtxs; i++ )
	{
		int j=hashtable[i];
		if ( j > -1 )
		{
			int k;
			for ( k=Mc->xadj[j]; k<Mc->xadj[j+1]; k++ )
			{
				Mr->adjncy[Mr_counter]=rmap[Mc->adjncy[k]];
				Mr->adjwgt[Mr_counter++]=Mc->adjwgt[k];
			}
		}
		Mr->xadj[i+1]=Mr_counter;
	}
	free(hashtable);

	printf("Projected %d flow values\n",Mr->nnz);
	printf("avg nnz per column:%d\n",(Mr_counter/Mr->nvtxs));

	return Mr;
}

/*
void mis_mlrmcl(GraphType* graph, idxtype* indices, Options opt)
{
	CtrlType c;
	int i,j,numLevels;
	Matrix** listOfMatrices;
	Matrix *M,*Mnew;
	idxtype* rmap=idxmalloc(graph->nvtxs, "mis_mlrmcl:rmap");

	c.CoarsenTo= opt.coarsenTo;

	InitRandom(time(NULL));
	
	listOfMatrices=mis_Coarsen2Way(&c,graph,opt.ncutify,&numLevels);

	printf("Num of levels of coarsening:%d\n",numLevels);
	
	M = listOfMatrices[numLevels-1];
	for( i=numLevels-1; i>0; i-- )
	{
		int nvtxs=listOfMatrices[i]->nvtxs;
		memcpy(rmap,listOfMatrices[i]->rmap,nvtxs*sizeof(idxtype));

		// listOfMatrices[i] will be freed in dprmcl. Hence
		// backup necessary info above. M will also be freed.
		M = dprmcl(M,listOfMatrices[i], graph, opt, 
					opt.iter_per_level,i);

		Mnew = mis_projectFlow(M,rmap,listOfMatrices[i-1]->nvtxs);
	
		if ( M != NULL )
		{
			freeMatrix(M);
		}
		M=Mnew;
	}

	M=dprmcl(M,listOfMatrices[0],graph,opt,opt.num_last_iter,0);

	getAttractorsForAll(M);
	idxcopy(M->nvtxs, M->attractors, indices);

	free(M);
	free(listOfMatrices);
	free(rmap);
}
*/

void mlmcl(int* nvtxs, idxtype* xadj, idxtype* adjncy, idxtype
*vwgt, idxtype* adjwgt, int* wgtflag, idxtype* indices, Options opt)
{
 /*	GraphType graph;
	my_SetUpGraph(&graph, *nvtxs, xadj, adjncy, vwgt, adjwgt,
	*wgtflag, 1); */
	int hubRemoval=opt.hubRemoval, recursiveCluster=0;
	float hub_pct = opt.hubPct;

	GraphType *graph = (GraphType*)malloc(sizeof(GraphType));
	my_SetUpGraph(graph, *nvtxs, xadj, adjncy, vwgt, adjwgt,
	*wgtflag, 1);
	// The last argument indicates we are setting up the original
	// graph 

	idxtype* newIds;
	if ( hubRemoval > 0 )
	{
		int hubThreshold = (int) floor(hub_pct * graph->nvtxs);
		GraphType *new_graph;
		newIds = removeHubs(graph, hubThreshold, *wgtflag,
						&new_graph, 0);
		free(graph->gdata);
		free(graph);
		graph = new_graph;
		
		// now need to remove any nodes that became singletons
		// because of hub removal.

		// we'll do another iteration of newIds, so back up 
		// the old newIds. newIds_bkp is of size *nvtxs.
		idxtype *newIds_bkp = newIds;

		int noOfSingletons = 0, newIdCounter;
		newIds = lookForSingletons(graph, &noOfSingletons);
		newIdCounter = graph->nvtxs - noOfSingletons;

		if ( noOfSingletons > 0 )
		{
			printf("%d nodes became singletons due to hub removal", 
						noOfSingletons );
			printf("; they will be removed.\n");
			fflush(stdout);

			getSubgraph(graph, newIds, newIdCounter, *wgtflag, 
							&new_graph);
			free(graph->gdata);
			free(graph);
			graph = new_graph;
			
			int i;
			for ( i=0; i<*nvtxs; i++ )
			{
				if ( newIds_bkp[i] > -1 )
				{
					newIds_bkp[i] = newIds[newIds_bkp[i]];
				}
				else
					newIds_bkp[i] = -1;
			}
			free(newIds);
		}

		newIds=newIds_bkp;
	}

//	printf("nnz:%d\n",graph.xadj[*nvtxs]);
	if ( opt.mis_coarsenType > 0 )
	{
//		mis_mlrmcl(graph, indices, opt); 
	}
	else
	{
		mlmclWithGraph(graph, indices, opt);
	}

	if ( hubRemoval > 0 )
	{
		int npart=mapPartition(indices, graph->nvtxs);
		float ncut=ComputeNCut(graph, indices, npart);
		printf("In graph that does not include hubs,"); 
		printf("No. of Clusters:%d, N-Cut value: %.2f\n", npart, ncut);

		mapIndices(indices, newIds, *nvtxs, npart);
		free(newIds);
		if ( *nvtxs - graph->nvtxs > 0 )
		{
			char filename[256];
			sprintf(filename, "input.nohubs.%.3f", hub_pct);
			WriteGraph(filename, graph->nvtxs, graph->xadj,
			graph->adjncy);
			printf("Wrote nohubs graph to %s\n", filename);
		}
	}

	if ( recursiveCluster > 0 )
	{
		int npart = mapPartition(indices, graph->nvtxs);
		float ncut = ComputeNCut(graph, indices, npart);
		printf("No. of clusters:%d, N-Cut:%.2f\n", npart, ncut);
		idxtype* hist = histogram(indices, graph->nvtxs, npart);

		int max=0, i=0, maxCluster=-1;
		for( i=0; i<npart; i++ )
		{
			if ( hist[i] > max )
			{
				max = hist[i];
				maxCluster = i;
			}
		}

		free(hist);

		if ( max > graph->nvtxs * 0.3 )
		{
			printf("Will recursively partition cluster of size");
			printf(" %d\n", max);
			
			idxtype* newIds = idxmalloc(graph->nvtxs,"mlmcl:newIds");
			int newIdCounter=0;
			for ( i=0; i<graph->nvtxs; i++ )
			{
				if ( indices[i] == maxCluster )
					newIds[i]=newIdCounter++;
				else
					newIds[i]=-1;
			}
			
			GraphType *new_graph;
			getSubgraph(graph, newIds, max, *wgtflag, &new_graph);

			idxtype *new_indices = idxmalloc(max,"mlmcl:new_indices");
			opt.coarsenTo = (int) round(((float) max 
							/ (float)graph->nvtxs) * opt.coarsenTo);
			mlmcl(&max,new_graph->xadj, new_graph->adjncy,
			new_graph->vwgt, new_graph->adjwgt, wgtflag,
			new_indices, opt );

			int new_npart = mapPartition( new_indices, max);
			for ( i=0; i<graph->nvtxs; i++ )
			{
				if ( newIds[i] > -1 )
				{
					int ni = new_indices[newIds[i]];
					if ( ni > 0 )
						indices[newIds[i]] = npart + ni - 1;
					else
						indices[newIds[i]] = maxCluster;
				}
			}
			
			printf("Recursive clustering yielded %d new",new_npart);
			printf(" clusters.");

			free(new_indices);
			free(newIds);
			free(new_graph->gdata);
			free(new_graph);

		}

	}
}


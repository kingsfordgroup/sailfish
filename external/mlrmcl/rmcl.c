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

void dprmclWrapper(int* nvtxs, idxtype* xadj, idxtype* adjncy, idxtype
*vwgt, idxtype* adjwgt, int* wgtflag, idxtype* indices, Options opt  )
{
	Matrix *flows, *adj, *final;
	GraphType graph;

	my_SetUpGraph(&graph, *nvtxs, xadj, adjncy, vwgt, adjwgt,
	*wgtflag, 1);
	
	adj=setupCanonicalMatrix(*nvtxs,graph.nedges,xadj,adjncy,adjwgt,opt.ncutify);
	flows=adj;

	final=dprmcl(flows,adj,&graph,opt,opt.num_last_iter,0);

	getAttractorsForAll(final);
	idxcopy(*nvtxs, final->attractors, indices);
	freeMatrix(final);

}

/*
 * Removes nodes which have locally converged from the graph.
 * i.e. if max(flows(:,i))>0.9, then adj(i,i)=1 and adj(j,i)=0,
 * for all j != i. 
 */
Matrix* changeAdjMatrix(Matrix* flows, Matrix* adj, idxtype*
			hashtable)
{
	int i,j,numConverged=0;
	Matrix* ret=allocMatrix(adj->nvtxs,adj->nnz,0,0,0);

	for (i=0; i<adj->nvtxs; i++)
	{
		hashtable[i]=-1;
	}

	ret->xadj[0]=0;
	for (i=0; i<adj->nvtxs; i++)
	{
		wgttype max=0;
		for (j=flows->xadj[i]; j<flows->xadj[i+1]; j++)
		{
			if ( flows->adjwgt[j] > max )
				max=flows->adjwgt[j];
		}

		if ( max > 0.9 )
		{
			ret->xadj[i+1]=ret->xadj[i]+1;
			ret->adjncy[ret->xadj[i]]=i;
			ret->adjwgt[ret->xadj[i]]=1;
			numConverged++;
		}
		else
		{
			int k=ret->xadj[i],l=adj->xadj[i];
			for( j=0; j<adj->xadj[i+1]-adj->xadj[i]; j++)
			{
				ret->adjncy[k+j]=adj->adjncy[l+j];
				ret->adjwgt[k+j]=adj->adjwgt[l+j];
			}
			ret->xadj[i+1]=ret->xadj[i]+j;
		}

	}

	freeMatrix(adj);
	printf("number of locally converged	nodes:%d\n",numConverged);
	return ret;
}

/* Warning: This function frees both the flows and the adj
 * matrices. */
Matrix* dprmcl(Matrix* flows, Matrix* adj, GraphType* graph, 
Options opt,  int maxiter, int level)
{
	Matrix *M0,*M, *Mnew;
	int numiter,numUnique,diff,i,numWithoutAttrs;
	idxtype* prevAttractors=idxmalloc(adj->nvtxs,"Attractors");
	idxtype* hashtable=idxmalloc(adj->nvtxs,"Hashtable");
	int calcChange=0,getNcut=0,dprPhase=0,calcDiff=0;
	float gamma = opt.gamma;
	int exact = opt.exact, k=opt.k;
	wgttype dpr_threshold = opt.dpr_threshold;
//	wgttype dpr_threshold = 1.5;
	wgttype changeInNorm=-1,changeInNcut=-1,prevNcut=-1;
	timer expTimer,infTimer,pnTimer,totalTimer;

	//printf("In dprmcl: penalty_power:%.1f\n", opt.penalty_power);
	// clear hashtable.
	for(i=0;i<adj->nvtxs;i++)
		hashtable[i]=-1;

	M0=adj;
	M=flows;

	/* Main Loop!! */
	numiter=0;
	Matrix* Mold=NULL;
	do
	{
		if ( dprPhase <= 0 && numiter > 9 && level <= 0 
			&& (changeInNorm != -1 && changeInNorm < 1.0e-03) )
		//	&& (changeInNcut != -1 && changeInNcut < 5.0e-03) )
		{
			dprPhase=1;
			//M0=getDprAdjMatrix(M,adj,hashtable,dpr_threshold);
//			printf("nnz of prior adj matrix:%d,",M0->nnz);
			M0=getDprAdjMatrix(M,M0,hashtable,dpr_threshold);
			adj=NULL;
//			printf("nnz of dpr adj matrix:%d\n",M0->nnz);
		}

/*		if ( dprPhase <= 0 && numiter > 4 )
		{
			M0=changeAdjMatrix(M,M0,hashtable);
			adj=NULL;
		} */

		if ( numiter > 9 && level <= 0 )
		{
			calcChange=1; //getNcut=1;
		}

		cleartimer(totalTimer);
		starttimer(totalTimer);

/*		if ( dprPhase > 0 )
		{
			M0=getDprAdjMatrix(M,adj,hashtable,dpr_threshold);
			printf("nnz of actual adj matrix:%d,",adj->nnz);
			printf("nnz of dpr adj matrix:%d\n",M0->nnz);
		}
*/		
		Mnew=allInOneStep(M,M0,hashtable,opt,calcChange,
		dprPhase);

/*		wgttype* rowSums = getRowSums(Mnew);
		int maxRowSumIndex = samax(Mnew->nvtxs, rowSums);
		wgttype maxRowSum = rowSums[maxRowSumIndex];
		printf("maxRowSum:%.3f, index of maxRowSum:%d\n",
				maxRowSum, maxRowSumIndex);
		free(rowSums);
*/
		if ( calcChange > 0 )
		{
			wgttype temp=changeInNorm;	
			changeBetweenMatrices(Mnew,M,&changeInNorm);///(M->nvtxs);
			changeInNorm=changeInNorm/M->nvtxs;
			if ( fabs(changeInNorm-temp) < 1.0e-04 )
			{
/*				char fname[256];
				sprintf(fname, "%d.M.indices", numiter);
				getAttractorsForAll(M);
				my_WritePartition(fname, M->attractors,
				M->nvtxs, 2);
				getAttractorsForAll(Mnew);
				sprintf(fname, "%d.Mnew.indices", numiter);
				my_WritePartition(fname, Mnew->attractors,
				Mnew->nvtxs, 2);
*/				
			/*	changeBetweenMatrices(Mnew, Mold, &temp);
				temp = temp/M->nvtxs;
				printf("Diff between Mnew & Mold:%.4f\n",
				temp);*/
			} 
		}

		stoptimer(totalTimer);

		if ( M != M0 )
		{
			freeMatrix(M);
			if ( flows != NULL )
				flows = NULL;
		}

/*		if ( Mold != M0 && Mold != NULL )
		{
			freeMatrix(Mold);
			if ( flows != NULL )
				flows = NULL;
		}
		Mold=M; */

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
			//	printf("diff:%d,numWithoutAttrs:%d\n",
			//		diff,numWithoutAttrs);
			}
			if ( calcChange > 0 )
			{}
		//		printf("avg. change:%.6f\n",changeInNorm);
		
/*			printf("Done with %d iterations\n", numiter);
			printf("Total time for iteration:%7.3f\n",
						gettimer(totalTimer)); */

			fflush(stdout);
		} 
		numiter++;

		if ( calcChange > 0 && changeInNorm < 1.0e-05 )
	//	if ( calcChange > 0 && changeInNorm < 1.0e-03 )
			break;

	}while(	numiter < maxiter );

	free(hashtable);
	free(prevAttractors);
	if ( adj != NULL )
	{
		if ( flows == adj )
			flows = NULL;
		freeMatrix(adj);
	}
	if ( M0 != NULL && M0 != adj )
	{
		freeMatrix(M0);
	}
	if ( flows != NULL )
	{
		freeMatrix(flows);
	}

	return M;
	
}

/*
void rmcl(int* nvtxs, idxtype* xadj, idxtype* adjncy, idxtype
*vwgt, idxtype* adjwgt, int* wgtflag, idxtype* indices,float
gamma, int ncutify,int maxiter, int exact, int k )
{
	Matrix *M0,*M, *Mnew;
	int nedges=xadj[*nvtxs], numiter,numUnique,diff,i;
	idxtype* prevAttractors=idxmalloc(*nvtxs,"Attractors");
	idxtype* hashtable=idxmalloc(*nvtxs,"Hashtable");
	int calcChange=0;

	// clear hashtable.
	for(i=0;i<*nvtxs;i++)
		hashtable[i]=-1;
	
	M0=allocMatrix(*nvtxs,nedges,1,0,0);
	M0->nvtxs=*nvtxs;
	M0->nnz=nedges;
	idxcopy(*nvtxs+1, xadj, M0->xadj); 
	idxcopy(nedges, adjncy, M0->adjncy);
	if ( (*wgtflag&1) == 1 )
	{
		int i;
		for(i=0;i<nedges;i++)
			M0->adjwgt[i]=(wgttype)adjwgt[i];
		normalizeColumns(M0,1,0);
	}
	else
	{
		if ( ncutify )
		{
			ncutifyWeights(M0,0,ncutify);
//			printf("adjwgt[30]:%1.4f\n",M0->adjwgt[30]);
			normalizeColumns(M0,1,0);
		}
		else
			normalizeColumns(M0,0,0);
	}

//	printf("adjwgt[30]:%1.4f\n",M0->adjwgt[30]);
	
	numiter=0;
	M=M0;
	do
	{
		if ( numiter > 0 )
		{
			idxcopy(*nvtxs,M->attractors,prevAttractors);
		}
			
//		Mnew=expandInflate(M,M0,gamma);
//		Mnew=expand(M,M0);
		Mnew=allInOneStep(M,M0,hashtable,gamma,exact,k,calcChange,0);
		if ( M != M0 )
		{
			freeMatrix(M);
		}
		M=Mnew;
	//	numUnique=noOfUnique(Mnew->attractors,Mnew->nvtxs,hashtable);

		numiter++;
		if ( numiter > 1 )
			diff=isConverged2(prevAttractors,M->attractors,M->nvtxs);
		else	
			diff=M->nvtxs;
	//	printf("no. of attractors:%d,diff:%d\n",numUnique,diff);

		if ( numiter%3 == 0 )
		{
			printf("Done with %d iterations, diff:%d\n",
			numiter,diff); 
			fflush(stdout);
		} 
//	}while( !isConverged(prevAttractors,M->attractors,*nvtxs) &&
//	numiter < maxiter );
	}while(diff>(0.001)*((double)M->nvtxs) && numiter<maxiter);

//	printf("Final diff:%d\n", diff);

	getAttractorsForAll(M);
	idxcopy(*nvtxs, M->attractors, indices);
	freeMatrix(M);
	Mnew=NULL; // M and Mnew are the same at this point.
	if ( M0 != NULL )
		freeMatrix(M0);
}
*/



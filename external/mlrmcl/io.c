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
 * io.c
 *
 * This file contains routines related to I/O
 *
 * Started 8/28/94
 * George
 *
 * $Id: io.c,v 1.1 2010-12-02 18:04:50 venu Exp $
 *
 */

#include <metis.h>

/************************************************************************
 * This function reads in a clustering file and return #clusters
 ***********************************************************************/
int readClustering(const char *filename, int *part, int nvtex){
  FILE *fpin;
  int i, numClus=0;

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("\n\nFailed to open clustering file %s\n", filename);
    exit(0);
  }

  for (i=0; i<nvtex; i++){
    fscanf(fpin, "%d", &part[i]);
    numClus = numClus > part[i]? numClus : part[i];
  }
  return numClus+1;
}

void ReadMatrix(Matrix *m, char *filename, wgttype threshold)
{
  idxtype *xadj, *adjncy;
  wgttype *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  sscanf(line, "%d %d", &(m->nvtxs), &(m->nnz));

  m->nnz *= 2;
  m->nnz++;
  m->xadj = idxmalloc(m->nvtxs, "ReadMatrix:xadj");
  m->adjncy = idxmalloc(m->nnz, "ReadMatrix:adjncy");
  m->adjwgt = (wgttype*)GKmalloc(sizeof(wgttype)*m->nnz, "ReadMatrix:adjwgt");

  long k=0, numPruned = 0;
  int edge,i;
  wgttype wt;
  for (m->xadj[0]=0, k=0, i=0 ; i<m->nvtxs; i++) 
  {
    
	do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));

	if ( feof(fpin) )
		break;
	
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

   
    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      wt = (wgttype)strtod(oldstr, &newstr);
      oldstr = newstr;

      if (edge < 0)
        break;

	  if ( wt < threshold )
	  {
	  	numPruned++;
		continue;
	  }

	  if ( k >= m->nnz )
	  {
	  	printf("No. of edges more than %d\n", m->nnz);
		printf("Reading list of vertex %d\n",i+1);
		abort();
	  }

      m->adjncy[k] = edge;
      m->adjwgt[k] = wt;
      k++;
    } 
    m->xadj[i+1] = k;
  }

  if ( threshold > 0 || numPruned > 0 )
  {
  	printf("Pruned %ld edges with wt less than %f\n", numPruned,
	threshold);
  }
  if ( k != m->nnz )
  {
  	printf("Promised edges: %d, Actual edges: %ld\n", m->nnz, k);
	m->adjncy = (idxtype*) realloc(m->adjncy, sizeof(idxtype)*k);
	m->adjwgt = (wgttype*) realloc(m->adjwgt, sizeof(wgttype)*k);
	m->nnz = k;
  }

  printf("Successfully read matrix\n");
  fflush(stdout);

  fclose(fpin);
  
}

void ReadGraph(GraphType *graph, const char *filename, int *wgtflag,
	int addSelfLoop, int txtFormat)
{
  int i, j, k, l, fmt, readew, readvw, ncon, edge, ewgt,
  actualnvtxs=0, noOfSelfLoops=0, allocedEdges; 
  
  idxtype *xadj, *adjncy, *vwgt, *adjwgt;
  char *line, *oldstr, *newstr;
  FILE *fpin;

  InitGraph(graph);

  line = (char *)malloc(sizeof(char)*(MAXLINE+1));

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  do {
    fgets(line, MAXLINE, fpin);
  } while (line[0] == '%' && !feof(fpin));

  if (feof(fpin)) {
    graph->nvtxs = 0;
    free(line);
    return;
  }

  fmt = ncon = 0;
  sscanf(line, "%d %d %d %d", &(graph->nvtxs), &(graph->nedges), &fmt, &ncon);
  
  readew = (fmt%10 > 0);
  readvw = ((fmt/10)%10 > 0);
  if (fmt >= 100) {
    printf("Cannot read this type of file format!");
    exit(0);
  }

  *wgtflag = 0;
  if (readew)
    *wgtflag += 1;
  if (readvw)
    *wgtflag += 2;

  if (ncon > 0 && !readvw) {
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("You specified ncon=%d, but the fmt parameter does not specify vertex weights\n", ncon);
    printf("Make sure that the fmt parameter is set to either 10 or 11.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  graph->nedges *=2;
  /* Venu: my addition */
  if ( addSelfLoop > 0 )
  	graph->nedges += graph->nvtxs;

  ncon = graph->ncon = (ncon == 0 ? 1 : ncon);

  /*printf("%d %d %d %d %d [%d %d]\n", fmt, fmt%10, (fmt/10)%10, ncon, graph->ncon, readew, readvw);*/

  if (graph->nvtxs > MAXIDX) 
    errexit("\nThe matrix is too big: %d [%d %d]\n", graph->nvtxs, MAXIDX, sizeof(idxtype));

  xadj = graph->xadj = idxsmalloc(graph->nvtxs+1, 0, "ReadGraph: xadj");
  adjncy = graph->adjncy = idxmalloc(graph->nedges, "ReadGraph: adjncy");

  vwgt = graph->vwgt = (readvw ? idxmalloc(ncon*graph->nvtxs, "ReadGraph: vwgt") : NULL);
  adjwgt = graph->adjwgt = (readew ? idxmalloc(graph->nedges, "ReadGraph: adjwgt") : NULL);

  allocedEdges = graph->nedges;
  noOfSelfLoops = 0;

  /* Start reading the graph file */
  for (xadj[0]=0, k=0, i=0, actualnvtxs=0; i<graph->nvtxs;i++) 
  {
	int vertexid;
    
	do {
      fgets(line, MAXLINE, fpin);
    } while (line[0] == '%' && !feof(fpin));

	if ( feof(fpin) )
		break;
	
    oldstr = line;
    newstr = NULL;

    if (strlen(line) == MAXLINE) 
      errexit("\nBuffer for fgets not big enough!\n");

    if (readvw) {
      for (l=0; l<ncon; l++) {
        vwgt[i*ncon+l] = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }
    }

	if ( txtFormat > 0 )
	{
		vertexid=(int)strtol(oldstr,&newstr,10)-1;
		oldstr=newstr;

		/* add missing nodes */
		while( i < vertexid )
		{
			if ( addSelfLoop > 0 )
			{
			/* self loop for missing nodes */
				xadj[i+1]=xadj[i]+1;
				adjncy[xadj[i]]=i;
				k++;
				if ( readew )
					adjwgt[xadj[i]]=1;
			}
			else
				xadj[i+1]=xadj[i];

			i++;
		}
		actualnvtxs++;
	}
	if ( addSelfLoop > 0 )
	{
		adjncy[k]=i; 
		/* For now, assigning an edge weight of 1 */
		// this weight is modified later on, being set to the
		// max. weight out of the node.
		if ( readew )
			adjwgt[k] = 1; 
		k++;
	}  

	idxtype maxwgt = 0;
	// keep track of max. edge weight to this node
    for (;;) {
      edge = (int)strtol(oldstr, &newstr, 10) -1;
      oldstr = newstr;

      if (readew) {
        ewgt = (int)strtol(oldstr, &newstr, 10);
        oldstr = newstr;
      }

      if (edge < 0)
        break;

	  /* Venu: my addition. Since we've already added self-loop
	   * above, disregard self loop here. */
	  if ( edge == i && addSelfLoop > 0 )
	  {
		graph->nedges--;
		noOfSelfLoops++;
		continue;
	  }
	  // if addSelfLoop <= 0, do not include self loops.
	  else if ( edge == i )
	  {
	  	graph->nedges--;
		noOfSelfLoops++;
		continue;
	  }

	  if ( k >= allocedEdges )
	  {
	  	printf("No. of edges more than %d\n", allocedEdges );
		printf("Reading list of vertex %d\n",i+1);
		abort();
	  }

      adjncy[k] = edge;
      if (readew) 
	  {
        adjwgt[k] = ewgt;
		if ( ewgt > maxwgt )
			maxwgt = ewgt;
	  }
      k++;
    } 
    xadj[i+1] = k;
	if ( addSelfLoop > 0 && readew )
	{
		// weight on the self-loop is set to the max. edge weight
		// out of the node.
		adjwgt[xadj[i]] = maxwgt;
	}
  }


  /* add missing nodes at the end of the file */
  while( i < graph->nvtxs )
  {
  	/* self loop for missing nodes */
  	if ( addSelfLoop > 0 )
  	{
  		xadj[i+1]=xadj[i]+1;
  		adjncy[xadj[i]]=i;
  		k++;
  		if ( readew )
  			adjwgt[xadj[i]]=1;
  	}
  	else
  		xadj[i+1]=xadj[i];
  	i++;
  }

  if ( addSelfLoop == 0  && noOfSelfLoops > 0 )
  	printf("removed %d self loops\n", noOfSelfLoops);

  if ( graph->nedges != k )
  {
	  printf("expected no. of edges:%d, xadj[nvtxs]:%d\n", 
  			graph->nedges, k);
  }
	

  if (k > allocedEdges) {
	  if ( addSelfLoop > 0 )
	  {
	  	k -= graph->nvtxs;
		graph->nedges -= graph->nvtxs;
	  }
	if ( txtFormat <= 0 )
	{
		k /= 2;
		graph->nedges /= 2;
	}
    printf("------------------------------------------------------------------------------\n");
    printf("***  I detected an error in your input file  ***\n\n");
    printf("In the first line of the file, you specified that the graph contained\n%d edges. However, I found %d edges in the file.\n", graph->nedges, k);
    if (2*k == graph->nedges) {
      printf("\n *> I detected that you specified twice the number of edges that you have in\n");
      printf("    the file. Remember that the number of edges specified in the first line\n");
      printf("    counts each edge between vertices v and u only once.\n\n");
    }
    printf("Please specify the correct number of edges in the first line of the file.\n");
    printf("------------------------------------------------------------------------------\n");
    exit(0);
  }

  if ( k != allocedEdges )
  {
//  	printf("Promised edges: %d, Actual edges: %ld\n", m->nnz, k);
	graph->adjncy = (idxtype*) realloc(graph->adjncy,sizeof(idxtype)*k);
	if ( readew )
	{
		graph->adjwgt = (idxtype*) realloc(graph->adjwgt,sizeof(idxtype)*k);
	}
  }

  graph->nedges = k;
  printf("Final number of edges is %d\n",graph->nedges);
  printf("First entry is %d\n",graph->xadj[0]);
  free(line);
}

void readMemberships( char *file, int numMemberships, idxtype
*members, idxtype **groups_xadj, idxtype **groupIds, idxtype
*ngroups)
{
	idxtype *groups = idxmalloc(numMemberships,
					"readMemberships:groups");
	
	FILE *fp; 
	if ( (fp = fopen(file, "r")) == NULL) 
	{
		printf("Failed to open file %s\n", file);
		exit(0);
	}

	int i,j;
	for ( i=0; i<numMemberships; i++ )
	{
		fscanf(fp, "%d %d", &members[i], &groups[i]);
		members[i]--; // should start at 0, not 1.
	}
	fclose(fp);

	ParallelQSortInts(groups, members, 0, numMemberships-1);

	*ngroups = 1;
	for ( i=1; i<numMemberships; i++ )
	{
		if ( groups[i] != groups[i-1] )
			(*ngroups)++;
	}

	*groups_xadj = idxmalloc( *ngroups + 1,
						"readMemberships:groups_xadj");
	*groupIds = idxmalloc( *ngroups, "readMemberships:groupIds");
	(*groups_xadj)[0] = 0;
	(*groupIds)[0] = groups[0];
	for ( i=1, j=0; i<numMemberships; i++ )
	{
		if ( groups[i] != groups[i-1] )
		{
			(*groupIds)[++j] = groups[i];
			(*groups_xadj)[j] = i;
		}
	}
	(*groups_xadj)[++j] = i;

	GKfree( (void**)&groups, LTERM);
}


/*************************************************************************
* This function writes out the partition vector, with gamma
* instead of nparts specified
**************************************************************************/
/* Venu: my addition*/
void my_WritePartition(const char *fname, idxtype *part, int n, float gamma)
{
  FILE *fpout;
  int i;

  if ((fpout = fopen(fname, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", fname);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);

}

void my_WritePartitionAddOne(const char *fname, idxtype *part, int n)
{
  FILE *fpout;
  int i;

  if ((fpout = fopen(fname, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", fname);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",(part[i]+1));

  fclose(fpout);

}

void my_WriteMappedPartition(const char *fname, idxtype *part, idxtype
*nodeMap, int n)
{
  FILE *fpout;
  int i;

  if ((fpout = fopen(fname, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", fname);

  for (i=0; i<n; i++)
  {
  	if ( nodeMap[i] > -1 )
	    fprintf(fpout,"%d %d\n",i+1,part[nodeMap[i]]);
  }

  fclose(fpout);

}


/* Assumes map is in increasing order,i.e. map[i] < map[j]
 * whenever i < j (except when map[j] == -1)
 */
void WriteRMap(const char *fname, idxtype *map, int n)
{
	int i,j;
	FILE *fp = fopen(fname, "w");
	for ( i=0,j=0; i<n; i++ )
	{
		while ( map[j] != i )
			j++;
		fprintf(fp, "%d %d\n", i+1, j+1);
	}
	fclose(fp);
}

/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePartition(const char *fname, idxtype *part, int n, int nparts)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.part.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n",part[i]);

  fclose(fpout);

}


/*************************************************************************
* This function writes out the partition vectors for a mesh
**************************************************************************/
void WriteMeshPartition(const char *fname, int nparts, int ne, idxtype *epart, int nn, idxtype *npart)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.epart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<ne; i++)
    fprintf(fpout,"%d\n", epart[i]);

  fclose(fpout);

  sprintf(filename,"%s.npart.%d",fname, nparts);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the partition file: %s", filename);

  for (i=0; i<nn; i++)
    fprintf(fpout,"%d\n", npart[i]);

  fclose(fpout);


}



/*************************************************************************
* This function writes out the partition vector
**************************************************************************/
void WritePermutation(const char *fname, idxtype *iperm, int n)
{
  FILE *fpout;
  int i;
  char filename[256];

  sprintf(filename,"%s.iperm",fname);

  if ((fpout = fopen(filename, "w")) == NULL) 
    errexit("Problems in opening the permutation file: %s", filename);

  for (i=0; i<n; i++)
    fprintf(fpout,"%d\n", iperm[i]);

  fclose(fpout);

}


int isGraphConnected(GraphType *graph)
{
	idxtype* visited = idxmalloc(graph->nvtxs, 
								"isGraphConnected:graph->nvtxs");
	int i;
	for ( i=0; i<graph->nvtxs; i++ )
		visited[i] = 0;

	int nVisited=0;

	dfTraversal(graph, 0, visited, &nVisited );

	if ( nVisited < graph->nvtxs )
	{
		printf("No. of vertices visited: %d\n", nVisited);
		return 0;
	}
	else
		return 1;
}

idxtype* getNodesToComponentMap(Matrix* graph1, int *numComponents,
wgttype minWgt)
{
	sortAdjLists(graph1->nvtxs, graph1->xadj, graph1->adjncy,
	graph1->adjwgt);
	Matrix *graph2, *graph;
	graph2 = getTranspose(graph1);
	graph = add(graph1, graph2);
	freeMatrix(graph2);
	idxtype *ret = idxmalloc(graph->nvtxs,
		"getNodesToComponentMap:ret");
	idxtype* visited = idxmalloc(graph->nvtxs, 
					"getNodesToComponentMap:visited");

	int i;
	for ( i=0; i<graph->nvtxs; i++ )
	{
		visited[i] = 0;
		ret[i] = -1;
	}
	
	int totalVisited = 0;
	*numComponents = 0;
	while ( totalVisited < graph->nvtxs )
	{
		for ( i=0; i<graph->nvtxs && visited[i] > 0 ; i++);
		if ( i == graph->nvtxs )
			break;

		int nVisited=0;
		dfTraversalMatrix(graph, i, visited, &nVisited, minWgt);
		totalVisited += nVisited;

		for ( i=0; i<graph->nvtxs; i++ )
		{
			// if node was visited in this round of dfs
			if ( ret[i] == -1 && visited[i] > 0 )
				ret[i] = *numComponents;
		}

		(*numComponents)++;
	}

	free( visited );
	freeMatrix(graph);

	return ret;
}

idxtype* compSizeDistribution(GraphType* graph, int
*numComponents)
{
	idxtype* cSizes = idxmalloc(graph->nvtxs,
						"compSizeDistribution:graph->nvtxs");	
	idxtype* visited = idxmalloc(graph->nvtxs, 
								"isGraphConnected:graph->nvtxs");
	
	int i;
	for ( i=0; i<graph->nvtxs; i++ )
	{
		visited[i] = 0;
		cSizes[i] = 0;
	}

	int totalVisited = 0;
	*numComponents = 0;
	while ( totalVisited < graph->nvtxs )
	{
		for ( i=0; i<graph->nvtxs && visited[i] > 0 ; i++);
		if ( i == graph->nvtxs )
			break;

		int nVisited;
		dfTraversal(graph, i, visited, &(cSizes[*numComponents]));
		totalVisited += cSizes[*numComponents];
		(*numComponents)++;
	}
	
	free(visited);
	cSizes = (idxtype*) realloc ( cSizes,
							(*numComponents)*sizeof(idxtype) );
	return cSizes;

}

/*************************************************************************
* This function checks if a graph is valid
**************************************************************************/
int CheckGraph(GraphType *graph)
{
  int i, j, k, l, nvtxs, err=0;
  idxtype *xadj, *adjncy, *adjwgt;

  nvtxs = graph->nvtxs;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  adjwgt = graph->adjwgt;


  for (i=0; i<nvtxs; i++) {
    for (j=xadj[i]; j<xadj[i+1]; j++) {
      k = adjncy[j];

      if (i == k) {
        printf("Vertex %d contains a self-loop (i.e., diagonal entry in the matrix)!\n", i);
        err++;
      }
      else {
        for (l=xadj[k]; l<xadj[k+1]; l++) {
          if (adjncy[l] == i) {
            if (adjwgt != NULL && adjwgt[l] != adjwgt[j]) {
              printf("Edges (%d %d) and (%d %d) do not have the same weight! %d %d\n", i,k,k,i, adjwgt[l], adjwgt[adjncy[j]]);
              err++;
            }
            break;
          }
        }
        if (l == xadj[k+1]) {
          printf("Missing edge: (%d %d)!\n", k, i);
          err++;
//		  abort();
        }
      }
    }
  }

  if (err > 0) 
    printf("A total of %d errors exist in the input file. Correct them, and run again!\n", err);

  return (err == 0 ? 1 : 0);
}


/*************************************************************************
* This function reads the element node array of a mesh
**************************************************************************/
idxtype *ReadMesh(const char *filename, int *ne, int *nn, int *etype)
{
  int i, j, k, esize;
  idxtype *elmnts;
  FILE *fpin;

  if ((fpin = fopen(filename, "r")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fscanf(fpin, "%d %d", ne, etype);

  switch (*etype) {
    case 1:
      esize = 3;
      break;
    case 2:
      esize = 4;
      break;
    case 3:
      esize = 8;
      break;
    case 4:
      esize = 4;
      break;
    default:
      errexit("Unknown mesh-element type: %d\n", *etype);
  }

  elmnts = idxmalloc(esize*(*ne), "ReadMesh: elmnts");

  for (j=esize*(*ne), i=0; i<j; i++) {
    fscanf(fpin, "%d", elmnts+i);
    elmnts[i]--;
  }

  fclose(fpin);

  *nn = elmnts[idxamax(j, elmnts)]+1;

  return elmnts;
}

void WriteMappedTxtGraphWithWts(const char *filename, int fullNvtxs, 
	idxtype *xadj, idxtype *adjncy, idxtype *adjwgt, 
	idxtype* map, int mappedNvtxs)
{
  int i, j;
  FILE *fpout;

  printf("Will write to file %s\n",filename);

  if ( adjwgt == NULL )
  {
  	printf("WriteGraphWithWts: adjwgt is null\n");
	abort();
  }

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 1\n", fullNvtxs, xadj[mappedNvtxs]/2);

  for (i=0; i<mappedNvtxs; i++) 
  {
	if ( xadj[i] == xadj[i+1] )
	{
		continue;
	}
	fprintf(fpout, "%d", (map[i]+1));
    for (j=xadj[i]; j<xadj[i+1]; j++)
	{
		int k=map[adjncy[j]];
    	fprintf(fpout, " %d %d", map[adjncy[j]]+1, adjwgt[j]);
	}
  	fprintf(fpout, "\n");
  }

  fclose(fpout);
}

void WriteMatrix(const char *filename, int nvtxs, idxtype *xadj,
idxtype *adjncy, wgttype *adjwgt)
{
  int i, j;
  FILE *fpout;

  printf("Will write to file %s\n",filename);

  if ( adjwgt == NULL )
  {
  	printf("WriteMatrix: adjwgt is null\n");
	abort();
  }

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 1", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
	if ( xadj[i] == xadj[i+1] )
	{
/*		printf("No neighbours for node %d\n", i);
		abort(); */
		continue;
	}
	fprintf(fpout, "%d %.6f", adjncy[xadj[i]]+1,adjwgt[xadj[i]]);
    for (j=xadj[i]+1; j<xadj[i+1]; j++)
      fprintf(fpout, " %d %.6f", adjncy[j]+1, adjwgt[j]);
  }
  fprintf(fpout, "\n");

  fclose(fpout);
}

/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteGraphWithWts(const char *filename, int nvtxs, idxtype *xadj,
idxtype *adjncy, idxtype *adjwgt)
{
  int i, j;
  FILE *fpout;

  printf("Will write to file %s\n",filename);

  if ( adjwgt == NULL )
  {
  	printf("WriteGraphWithWts: adjwgt is null\n");
	abort();
  }

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 1", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
	if ( xadj[i] == xadj[i+1] )
	{
/*		printf("No neighbours for node %d\n", i);
		abort(); */
		continue;
	}
	fprintf(fpout, "%d %d", adjncy[xadj[i]]+1,adjwgt[xadj[i]]);
    for (j=xadj[i]+1; j<xadj[i+1]; j++)
      fprintf(fpout, " %d %d", adjncy[j]+1, adjwgt[j]);
  }
  fprintf(fpout, "\n");

  fclose(fpout);
}

/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteTxtGraph(const char *filename, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d", nvtxs, xadj[nvtxs]);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n%d", i+1);
    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %d", adjncy[j]+1);
  }
  fprintf(fpout, "\n");

  fclose(fpout);
}

/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteGraph(const char *filename, int nvtxs, idxtype *xadj, idxtype *adjncy)
{
  int i, j;
  FILE *fpout;

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d\n", nvtxs, xadj[nvtxs]/2);

  for (i=0; i<nvtxs; i++) {
  	if ( xadj[i] < xadj[i+1] )
	{
	  if ( adjncy[xadj[i]] != i )
		  fprintf(fpout, "%d", adjncy[xadj[i]]+1);
   	  for (j=xadj[i]+1; j<xadj[i+1]; j++)
	  {
       if ( adjncy[j] != i )
	      fprintf(fpout, " %d", adjncy[j]+1);
	  }
	}

    fprintf(fpout, "\n");
  }

  fclose(fpout);
}


/*************************************************************************
* This function writes a graphs into a file 
**************************************************************************/
void WriteMocGraph(GraphType *graph)
{
  int i, j, nvtxs, ncon;
  idxtype *xadj, *adjncy;
  float *nvwgt;
  char filename[256];
  FILE *fpout;

  nvtxs = graph->nvtxs;
  ncon = graph->ncon;
  xadj = graph->xadj;
  adjncy = graph->adjncy;
  nvwgt = graph->nvwgt;

  sprintf(filename, "moc.graph.%d.%d", nvtxs, ncon);

  if ((fpout = fopen(filename, "w")) == NULL) {
    printf("Failed to open file %s\n", filename);
    exit(0);
  }

  fprintf(fpout, "%d %d 10 1 %d", nvtxs, xadj[nvtxs]/2, ncon);

  for (i=0; i<nvtxs; i++) {
    fprintf(fpout, "\n");
    for (j=0; j<ncon; j++)
      fprintf(fpout, "%d ", (int)((float)10e6*nvwgt[i*ncon+j]));

    for (j=xadj[i]; j<xadj[i+1]; j++)
      fprintf(fpout, " %d", adjncy[j]+1);
  }

  fclose(fpout);
}

void printHistogram(idxtype* hist, int maxDegree, FILE* fp)
{
	int i;
	for ( i=0; i<maxDegree+1; i++ )
	{
		if ( hist[i] > 0 )
			fprintf(fp, "%d %d\n", i, hist[i]);
	}
}

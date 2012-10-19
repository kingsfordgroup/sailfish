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
 * proto.h
 *
 * This file contains header files
 *
 * Started 10/19/95
 * George
 *
 * $Id: proto.h,v 1.38 2010-11-30 23:19:00 venu Exp $
 *
 */

/* bucketsort.c */
void BucketSortKeysInc(int, int, idxtype *, idxtype *, idxtype *);

/* ccgraph.c */
void my_CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraph(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraphNoMask(CtrlType *, GraphType *, int, idxtype *, idxtype *);
void CreateCoarseGraph_NVW(CtrlType *, GraphType *, int, idxtype *, idxtype *);
GraphType *my_SetUpCoarseGraph(GraphType *, int);
GraphType *SetUpCoarseGraph(GraphType *, int, int);
void ReAdjustMemory(GraphType *, GraphType *, int);
void CreateCoarseGraph_PowerLaw(CtrlType*, GraphType*, int,
		idxtype*, idxtype*, idxtype* );

/* Venu: My addition */
/* mis_coarsen.c */
Matrix** mis_Coarsen2Way(CtrlType *, GraphType*, int, int*);

/* coarsen.c */
GraphType *Coarsen2Way(CtrlType *, GraphType *);

/* metrics.c */
double ahkEdges(GraphType*, float);
int numberOfMatchingEdges(GraphType*, GraphType*);
float stdDeviation(idxtype*, int);
idxtype* getWeightsHistogram(GraphType*, int*, int);
idxtype* getDegreeHistogram(GraphType*, int*, int);
void ComputeAdjWgtSums(GraphType*);
int ComputeNumEdgesCut(GraphType*, idxtype*);
idxtype ComputeVolSubset(GraphType*, idxtype*);
int ComputeCut(GraphType *, idxtype *);
float ComputeConductance(GraphType *, idxtype *, int );
float ComputeNCut(GraphType *, idxtype *, int );
float Overlap_ComputeNcutVector(GraphType*, idxtype*, idxtype*,
idxtype, float*);
float ComputeNCutVector(GraphType*, idxtype*, int, float*);
int mapPartition(idxtype*, idxtype);
idxtype* histogram(idxtype*, int, int);
//long* histogramForLongs(long*, long, long);
int checkValidUndirectedGraph(GraphType*);
int checkEqualityOfGraphs(GraphType*);

/* graph.c */
void my_SetUpGraph(GraphType*, int, idxtype*, idxtype*, idxtype*,
idxtype*, int, int);

/* match.c */
void Match_RM(CtrlType *, GraphType *);
void Match_RM_NVW(CtrlType *, GraphType *);
void Match_HEM(CtrlType *, GraphType *);
void Match_SHEM(CtrlType *, GraphType *);
void my_Match_HEMN(CtrlType *, GraphType *);
void my_Match_SHEMN(CtrlType *, GraphType *);
void my_Match_RM(CtrlType *, GraphType *);
void my_Match_PowerLaw_FC(CtrlType *, GraphType *);


/* mclutils.c */
void removeSelfLoops(GraphType*);
Matrix* removeSelfLoops(Matrix*, int);
void testCoarsening(int*, idxtype*, idxtype*, idxtype*, idxtype*,
int*, int);
void dump_graph(GraphType*);
void ncutifyWeights(Matrix*, int, int);
void normalizeColumns(Matrix*, int, int);
void printRow(Matrix*, int);
int numWithoutAttractors(Matrix*);
void getAttractorsForAll(Matrix*);
int isConverged(idxtype*, idxtype*, int);
int isConverged2(idxtype*, idxtype*, int);
int noOfUnique(idxtype*, int, idxtype*);
void freeMatrix(Matrix*);
Matrix* allocMatrix(int, int, int, int, int);
void dumpMatrix(Matrix*);
Matrix* setupCanonicalMatrix(int, int, idxtype*, idxtype*,
idxtype*, int );
void sortAdjLists(int, idxtype*, idxtype*, wgttype*);
void getPermutedGraph(idxtype*, idxtype*, int, int, idxtype*,
	idxtype*, idxtype*, idxtype**, idxtype**, idxtype** );
Matrix* permuteRowsAndColumns(Matrix*, idxtype*, idxtype*);

/* rmcl.c */
/*void rmcl(int*, idxtype*, idxtype*, idxtype*, idxtype*, 
int*, idxtype*,float,int,int, int, int);
*/
Matrix* dprmcl(Matrix*, Matrix*, GraphType*, Options, int, int);
void dprmclWrapper(int*, idxtype*,idxtype*,idxtype*, 
idxtype*, int*, idxtype*, Options);
 

/* mlmcl.c */
void mlmcl(int*, idxtype*, idxtype*, idxtype*, idxtype*, int*,
				idxtype*, Options);

/* pagerank.c */
void initPageRankOptions(PageRankOptions *opt);
wgttype* pagerank(GraphType*, PageRankOptions, Matrix *M = NULL); 


/* mclbase.c */
void exactPruneGraphCountingSort(GraphType*, int, int);
int sizeOfSetIntersect(idxtype*, int, idxtype*, int);
void knnMatrix(Matrix*, int);
void checksums(Matrix*);
Matrix* getTranspose(Matrix*);
Matrix* getTranspose2(Matrix*);
Matrix* AtimesATransposeBlocking(Matrix* , wgttype*, wgttype*,
					int, wgttype);
Matrix* AtimesATranspose(Matrix*, wgttype*, wgttype*, int, int,
			int, int, wgttype, timer*);
Matrix* expand(Matrix*,Matrix*);
Matrix* expand_ht(Matrix*,Matrix*,idxtype*, int, wgttype
threshold = 0);
Matrix* allInOneStep(Matrix*, Matrix*, idxtype*,Options ,int, int);
Matrix* getDprAdjMatrix(Matrix*,Matrix*,idxtype*,wgttype);
Matrix* add(Matrix*, Matrix*);
void changeBetweenMatrices(Matrix*,Matrix*,wgttype*);
void changeBetweenMatrices2(Matrix*,Matrix*,wgttype*);
void inflate(Matrix*,float);
void pruneAndNormalize(Matrix*,int,int);
wgttype exactPrune(int*, idxtype*, wgttype*, int, wgttype);
void exactPruneMatrix(Matrix*, int);
void exactPruneGraph(GraphType*, int, idxtype thresh=0);
void retainTopNeighborsPerNode(GraphType*, Hashtable* , float );
int compareints(const void *, const void *);

/* memory.c */
void my_AllocateWorkSpace(CtrlType *ctrl, GraphType *graph);
void AllocateWorkSpace(CtrlType *, GraphType *, int);
void FreeWorkSpace(CtrlType *, GraphType *);
int WspaceAvail(CtrlType *);
idxtype *idxwspacemalloc(CtrlType *, int);
void idxwspacefree(CtrlType *, int);
float *fwspacemalloc(CtrlType *, int);
void fwspacefree(CtrlType *, int);
GraphType *CreateGraph(void);
void InitGraph(GraphType *);
void FreeGraph(GraphType *);

/* myqsort.c */
void iidxsort(int, idxtype *);
void iintsort(int, int *);
void ikeysort(int, KeyValueType *);
void ikeyvalsort(int, KeyValueType *);

/* subgraphs.c */
idxtype* removeHubs(GraphType*, int, int, GraphType**, int);
void getSubgraph(GraphType*, idxtype*, int, int, GraphType**);
GraphType* getCutGraph(GraphType*, idxtype*, int);
void globallySampleEdges(int, int, idxtype*, idxtype*, idxtype*,
		idxtype**, idxtype**, float);

/* timing.c */
void InitTimers(CtrlType *);
void PrintTimers(CtrlType *);
double seconds(void);

/* util.c */
void assignClustersToHubs(idxtype*, idxtype*, int, int,
GraphType*);
void dfTraversalMatrix(Matrix*, idxtype, idxtype*, int*, wgttype
minWgt);
void dfTraversal(GraphType*, idxtype, idxtype*, int*);
void initOptions(Options*);
void mapIndices(idxtype*, idxtype*, int, int);
idxtype* lookForSingletons(GraphType*, int*);
void errexit(char *,...);
#ifndef DMALLOC
int *imalloc(int, const char *);
idxtype *idxmalloc(int, const char *);
float *fmalloc(int, const char *);
int *ismalloc(int, int, const char *);
idxtype *idxsmalloc(int, idxtype, const char *);
//void *GKmalloc(int, char *);
void *GKmalloc(long, const char *);
#endif
void GKfree(void **,...); 
int *iset(int n, int val, int *x);
idxtype *idxset(int n, idxtype val, idxtype *x);
float *sset(int n, float val, float *x);
int iamax(int, int *);
int idxamax(int, idxtype *);
int idxamax_strd(int, idxtype *, int);
int samax(int, float *);
int samax2(int, float *);
int idxamin(int, idxtype *);
int samin(int, float *);
int idxsum(int, idxtype *);
int idxsum_strd(int, idxtype *, int);
void idxadd(int, idxtype *, idxtype *);
int charsum(int, const char *);
int isum(int, int *);
float ssum(int, float *);
float ssum_strd(int n, float *x, int);
void sscale(int n, float, float *x);
float snorm2(int, float *);
float sdot(int n, float *, float *);
void saxpy(int, float, float *, int, float *, int);
void ParallelQSortFloatsInts(wgttype*, idxtype*, int, int);
void ParallelQSortIntsUsingScores(idxtype*, idxtype*, idxtype*,
			int, int);
void ParallelQSort(idxtype*,wgttype*,int,int);
void ParallelQSortInts(idxtype*,idxtype*,int,int);
void QSortIntsUsingInts(idxtype*, idxtype*, int, int);
void ParallelQSortLongs(long*,wgttype*,int,int);
int bsearch_insertPos(idxtype*, int, int, int);
void RandomPermute(int, idxtype *, int);
void permuteDegreeOrder(int, idxtype*, idxtype*);
wgttype RandomSelect(wgttype*, int, int, int);
idxtype RandomSelectInts(idxtype*, int, int, int);
double drand48();
void srand48(long);
int ispow2(int);
void InitRandom(int);
int log2(int);

/* io.c */
void ReadMatrix(Matrix *, char*, wgttype threshold=0);
void readMemberships(char*, int, idxtype*, idxtype**, idxtype**,
				idxtype*);
idxtype* getNodesToComponentMap(Matrix*, int*, wgttype);
idxtype* compSizeDistribution(GraphType*, int*);
int isGraphConnected(GraphType*);
void WriteRMap(const char*, idxtype*, int);
void printHistogram(idxtype*, int, FILE*);
int readClustering(const char *, int *, int);
void ReadGraph(GraphType *, const char *, int *, int, int);
//void ReadTxtGraph(GraphType *, char *, int *, int);
void WritePartition(const char *, idxtype *, int, int);
void my_WritePartition(const char *, idxtype *, int, float);
void my_WritePartitionAddOne(const char *, idxtype *, int);
void my_WriteMappedPartition(const char *, idxtype *, idxtype *, int);
void WriteMeshPartition(const char *, int, int, idxtype *, int, idxtype *);
void WritePermutation(const char *, idxtype *, int);
int CheckGraph(GraphType *);
idxtype *ReadMesh(const char *, int *, int *, int *);
void WriteGraph(const char *, int, idxtype *, idxtype *);
void WriteTxtGraph(const char *, int, idxtype *, idxtype *);
void WriteMatrix(const char *, int, idxtype*, idxtype*, wgttype*);
void WriteGraphWithWts(const char *, int, idxtype*, idxtype*,
idxtype*);
void WriteMappedTxtGraphWithWts(const char *, int,
		idxtype*,idxtype*,idxtype*,idxtype*, int);


/* merge.c */
ListGraph* createClusterGraph(const idxtype*, int, const GraphType* );
void mergeBestClusters(ListGraph*, idxtype*, int, int);

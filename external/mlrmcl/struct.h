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
 * struct.h
 *
 * This file contains data structures for ILU routines.
 *
 * Started 9/26/95
 * George
 *
 * $Id: struct.h,v 1.9 2010-11-30 23:18:37 venu Exp $
 */

/* Undefine the following #define in order to use short int as the idxtype */
#define IDXTYPE_INT

/* Indexes are as long as integers for now */
#ifdef IDXTYPE_INT
typedef int idxtype;
#else
typedef short idxtype;
#endif

/* Venu : my addition below */
/* Undefine WGTTYPE_FLOAT in order to use double as the wgttype
 * */
#define WGTTYPE_FLOAT
/* wgttype is used to store weights in the matrices*/
#ifdef WGTTYPE_FLOAT
typedef float wgttype;
#else
typedef double wgttype;
#endif

#define MAXIDX	(1<<8*sizeof(idxtype)-2)


/*************************************************************************
* The following data structure stores key-value pair
**************************************************************************/
struct KeyValueType {
  idxtype key;
  idxtype val;
};

typedef struct KeyValueType KeyValueType;


/*************************************************************************
* The following data structure will hold a node of a doubly-linked list.
**************************************************************************/
struct ListNodeType {
  int id;                       	/* The id value of the node */
  struct ListNodeType *prev, *next;     /* It's a doubly-linked list */
};

typedef struct ListNodeType ListNodeType;



/*************************************************************************
* The following data structure is used to store the buckets for the 
* refinment algorithms
**************************************************************************/
struct PQueueType {
  int type;                     /* The type of the representation used */
  int nnodes;
  int maxnodes;
  int mustfree;

  /* Linear array version of the data structures */
  int pgainspan, ngainspan;     /* plus and negative gain span */
  int maxgain;
  ListNodeType *nodes;
  ListNodeType **buckets;

  /* Heap version of the data structure */
  KeyValueType *heap;
  idxtype *locator;
};

typedef struct PQueueType PQueueType;


/*************************************************************************
* The following data structure stores an edge
**************************************************************************/
struct edegreedef {
  idxtype pid;
  idxtype ed;
};
typedef struct edegreedef EDegreeType;


/*************************************************************************
* The following data structure stores an edge for vol
**************************************************************************/
struct vedegreedef {
  idxtype pid;
  idxtype ed, ned;
  idxtype gv;
};
typedef struct vedegreedef VEDegreeType;


/*************************************************************************
* This data structure holds various working space data
**************************************************************************/
struct workspacedef {
  idxtype *core;			/* Where pairs, indices, and degrees are coming from */
  int maxcore, ccore;

  EDegreeType *edegrees;
  VEDegreeType *vedegrees;
  int cdegree;

  idxtype *auxcore;			/* This points to the memory of the edegrees */

  idxtype *pmat;			/* An array of k^2 used for eliminating domain 
                                           connectivity in k-way refinement */
};

typedef struct workspacedef WorkSpaceType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct rinfodef {
 int id, ed;            	/* ID/ED of nodes */
 int ndegrees;          	/* The number of different ext-degrees */
 EDegreeType *edegrees;     	/* List of edges */
};

typedef struct rinfodef RInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* vol-based partition
**************************************************************************/
struct vrinfodef {
 int id, ed, nid;            	/* ID/ED of nodes */
 int gv;            		/* IV/EV of nodes */
 int ndegrees;          	/* The number of different ext-degrees */
 VEDegreeType *edegrees;     	/* List of edges */
};

typedef struct vrinfodef VRInfoType;


/*************************************************************************
* The following data structure holds information on degrees for k-way
* partition
**************************************************************************/
struct nrinfodef {
 idxtype edegrees[2];  
};

typedef struct nrinfodef NRInfoType;

/* Venu: my addition */
/****************
* This data structures holds a (square) sparse matrix
*/
struct matrixdef{
	int nvtxs, nnz; /* nnz stands for number of non-zero
	entries*/
	idxtype* xadj; /* xadj[i+1]-xadj[i] gives the number of
	non-zero entries in column i */
	idxtype* adjncy;
	wgttype* adjwgt; /* array that stores weights of the
	adjacency lists */
	wgttype* adjwgtsum; /* sum of adjacency weights of each node,
	or the sum of a column, basically. assigned only during the
	expandAndInflate phase */
	wgttype* maxwgt; /* max wgt of each column, assigned only
	during the expandAndInflate phase */
	idxtype* attractors; /* the row with the maximum weight in a
	column, and which has weight > 0.50 */

	idxtype* rmap; /* This is for mis_coarsen */
	int currentSize;
	int sizeIncrement;
};

typedef struct matrixdef Matrix;

struct hashtabledef
{
	int numHashes;
	int numNodes;

	// sortedNodeIds is the list of nodeIds sorted according to
	// their fingerprints.
	idxtype *sortedNodeIds; 

	idxtype *hashes;
};

typedef struct hashtabledef Hashtable;

typedef unsigned int UInt32;
//typedef int UInt32;

struct threadData {
	int nvtxs;
	idxtype *xadj, *adjncy, *adjwgt;
	float threshold;
	idxtype **ret_xadj, **ret_adjncy, **ret_adjwgt;
	Matrix **ret;
};

/*************************************************************************
* This data structure holds the input graph
**************************************************************************/
struct graphdef {
  idxtype *gdata, *rdata;	/* Memory pools for graph and refinement data.
                                   This is where memory is allocated and used
                                   the rest of the fields in this structure */

  int nvtxs, nedges;		/* The # of vertices and edges in the graph */
  idxtype *xadj;		/* Pointers to the locally stored vertices */
  idxtype *vwgt;		/* Vertex weights */
  idxtype *vsize;		/* Vertex sizes for min-volume formulation */
  idxtype *adjncy;		/* Array that stores the adjacency lists of nvtxs */
  idxtype *adjwgt;		/* Array that stores the weights of the adjacency lists */

  idxtype *adjwgtsum;		/* The sum of the adjacency weight of each vertex */

  idxtype *label;

  idxtype *cmap;
  
  /* Venu: my addition
  * Refine map: maps vertices in coarse graph to vertices in
   * refined graph. Two maps needed as 2 vertices are mapped to
   * one vertex in the coarse graph */
  idxtype *rmap1;
  idxtype *rmap2;
  idxtype *numDescendants;
  /* indicates if this graph is the original
  (i.e. the most refined) graph*/
  int isOrgGraph; 
  int isDirected; 
  wgttype *pagerank;

  /* Partition parameters */
  int mincut, minvol;
  idxtype *where, *pwgts;
  int nbnd;
  idxtype *bndptr, *bndind;

  /* Bisection refinement parameters */
  idxtype *id, *ed;

  /* K-way refinement parameters */
  RInfoType *rinfo;

  /* K-way volume refinement parameters */
  VRInfoType *vrinfo;

  /* Node refinement information */
  NRInfoType *nrinfo;


  /* Additional info needed by the MOC routines */
  int ncon;			/* The # of constrains */ 
  float *nvwgt;			/* Normalized vertex weights */
  float *npwgts;		/* The normalized partition weights */

  struct graphdef *coarser, *finer;
};

typedef struct graphdef GraphType;

struct intlistdef{
	idxtype* l;
	int length;
	int allocSize;
	int increment;
};
typedef struct intlistdef ListInt;

struct wgtlistdef{
	wgttype* l;
	int length;
	int allocSize;
	int increment;
};
typedef struct wgtlistdef ListWgt;

struct listgraphdef{
	int nvtxs;
	int nedges;
	ListInt* adjLists;
	ListWgt* wgtLists;
	wgttype* vols; // volumes, i.e. sum of edge weights
};
typedef struct listgraphdef ListGraph;

/*************************************************************************
* The following data type implements a timer
**************************************************************************/
typedef double timer;


/*************************************************************************
* The following structure stores information used by Metis
**************************************************************************/
struct controldef {
  int CoarsenTo;		/* The # of vertices in the coarsest graph */
  int dbglvl;			/* Controls the debuging output of the program */
  int CType;			/* The type of coarsening */
  int IType;			/* The type of initial partitioning */
  int RType;			/* The type of refinement */
  int maxvwgt;			/* The maximum allowed weight for a vertex */
  float nmaxvwgt;		/* The maximum allowed weight for a vertex for each constrain */
  int optype;			/* Type of operation */
  int pfactor;			/* .1*prunning factor */
  int nseps;			/* The number of separators to be found during multiple bisections */
  int oflags;

  WorkSpaceType wspace;		/* Work Space Informations */

  /* Various Timers */
  timer TotalTmr, InitPartTmr, MatchTmr, ContractTmr, CoarsenTmr, UncoarsenTmr, 
        SepTmr, RefTmr, ProjectTmr, SplitTmr, AuxTmr1, AuxTmr2, AuxTmr3, AuxTmr4, AuxTmr5, AuxTmr6;

};

typedef struct controldef CtrlType;

struct optionsdef{
	int coarsenTo;
	float gamma;
	int iter_per_level;
	int num_last_iter;
	int ncutify;
	int hubRemoval;
	float hubPct;
	int mis_coarsenType;
	int matchType;
	int exact;
	int k;
	float dpr_threshold;
	int transformAdj;
	float penalty_power; // used in transformAdj
};

typedef struct optionsdef Options;

struct pagerankOptionsDef{
	wgttype alpha; // random jump probability.
	int max_iters;
	wgttype convergeThreshold;
};

typedef struct pagerankOptionsDef PageRankOptions;

struct dirToUndirOptionsDef{
	int conversionMethod;
	wgttype threshold;
	int blockSize;
	wgttype bcWeight;
	wgttype prunePercent;
	wgttype scale;
	int invDegreeType;
};

typedef struct dirToUndirOptionsDef DirToUndirOptions;

/*************************************************************************
* The following data structure stores max-partition weight info for 
* Vertical MOC k-way refinement
**************************************************************************/
struct vpwgtdef {
  float max[2][MAXNCON];
  int imax[2][MAXNCON];
};

typedef struct vpwgtdef VPInfoType;





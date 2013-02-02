

/*
 * This program is free software distributed under the LGPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/lgpl-3.0.txt
 */

/* tsnnls_test is a general command-line utility providing access to
   some of the functionality of the tsnnls library. It is mostly useful
   for debugging tsnnls, but it can be used to quickly check the results
   of a MATLAB or Octave computation from the command line. */

/*
 * Note that for both formats, row and column indexing is 1-based.
 *
 * Sparse format is a plain text file containing:
 *
 * SPARSE
 * row_number column_number
 * number_of_nonzeros
 * r1 c1 val11
 * r1 c2 val12
 * ...
 * r2 c1 val21
 * ...
 *
 * Dense format is a plain text file corresponding to the Octave text
 * format:
 *
 * Any number of lines starting with #, as long as they contain
 *
 * # rows: m
 * # columns: n
 *
 * val1_1 val1_2 val1_3 val1_4 ... val1_col_number
 * val2_1 ...
 * ...
 * val_row_number_1 ...
 *
 * where valm_n is the value of the mth row and nth column.
 *
 *
 * Representative files for both of these formats are included in the
 * directory 'test_files' in the distribution root directory.
 */

#include <config.h>

#define PACKAGE_VERSION "2.3.3"

#ifdef HAVE_CTYPE_H
  #include <ctype.h>
#endif

#include "lsqr.h"
#include "tsnnls.h"
#include "float.h"
//#ifdef HAVE_CBLAS_H
//  #include <cblas.h>
//#else
//  #ifdef HAVE_VECLIB_CBLAS_H
//    #include <vecLib/cblas.h>
//  #else
//    #ifdef HAVE_ATLAS_CBLAS_H
//      #include <atlas/cblas.h>
//    #endif
//  #endif
//#endif

//#ifdef HAVE_DARWIN          /* We use the Apple BLAS/LAPACK if possible. */
// #include <vecLib/vBLAS.h>
//#else
// #include "tsnnls/gsl_cblas.h"
//#endif

#ifdef HAVE_CLAPACK_H
  #include "clapack.h"
#else
  #ifdef HAVE_VECLIB_CLAPACK_H
    #include "vecLib/clapack.h"
  #else
    #ifdef HAVE_ATLAS_CLAPACK_H
      #include "atlas/clapack.h"
    #endif
  #endif
#endif

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef HAVE_SYS_TIME_H
  #include <sys/time.h>
#endif

#ifdef HAVE_SYS_RESOURCE_H
  #include <sys/resource.h>
#endif

#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#include "tsnnls_blas_wrappers.h"

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#include <argtable2.h>

int VERBOSITY = 0;
enum STYPE { tsnnls, pjv, tlsqr, fallback, SOLlsqr, spiv } solvers[10] = { tsnnls };
int nsolvers = 0;
int nconstrained = -1;
int exit_status = 0;

int read_sparse( FILE *fp, double **vals, int *dim, int *cols );
taucs_ccs_matrix* read_sparse_non_stupid( const char* fname, int flags);
int read_mat( FILE *fp, double **vals, int *dim, int *cols);

#ifdef HAVE_SYS_RESOURCE_H
struct rusage start, end;
#endif

int tpStorage;

static double* lsqrwrapper( taucs_ccs_matrix* Af, double* b );

static void
start_clock()
{
#ifdef HAVE_SYS_RESOURCE_H
	getrusage(RUSAGE_SELF, &start);
#endif
}

static double
end_clock()
{
#ifdef HAVE_SYS_RESOURCE_H
    getrusage(RUSAGE_SELF, &end);

    double diff = end.ru_utime.tv_sec+end.ru_utime.tv_usec*1e-6 -
                  start.ru_utime.tv_sec+start.ru_utime.tv_usec*1e-6;
	return diff;
#else
	// This timer doesn't currently work without rusage
	return -1;
#endif
}

double *loadvals(const char *filename,int *dim,int *cols)

/* Attempts to load A from filename. */

{

  FILE *Af;
  double *vals;

  Af = fopen(filename,"r");

  if (Af == NULL) {

    fprintf(stderr,"tsnnls_test: Couldn't open file %s.\n", filename);
    exit(1);

  }

  if (read_mat(Af, &vals, dim, cols) == -1) { /* If reading as dense fails. */

    if (read_sparse(Af, &vals, dim, cols)) { /* If read as sparse fails. */

      fprintf(stderr,"tsnnls_test: Couldn't parse file %s.\n", filename);
      exit(1);

    }

  }

  fclose(Af);

  return vals;

}

taucs_ccs_matrix*
read_sparse_non_stupid( const char* fname, int flags) {

  FILE* fp = fopen(fname,"r");
  int i,j;
  float val;
  int nnz;
  int nr,nc;
  int ctr;
  int read;

  read = fscanf(fp, "%d %d\n", &nr, &nc);
  read = fscanf(fp, "%d\n", &nnz);

  taucs_ccs_matrix* mat;// = taucs_ccs_create(nr, nc, nnz, flags);

  mat = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  mat->n = nc;
  mat->m = nr;
  mat->flags = flags;

  mat->colptr = (int*)malloc(sizeof(int)*(mat->n+1));
  mat->rowind = (int*)malloc(sizeof(int)*nnz);
  mat->values.s = (float*)malloc(sizeof(float)*nnz);
    /*
  vals = (double*)calloc(nnz,sizeof(double));
  colels = (double**)calloc(result->n,sizeof(double*));
  colrows = (int**)calloc(result->n,sizeof(int*));
  colcounts = (int*)calloc(result->n,sizeof(int));
    */

  printf("New (%d x %d) Matrix with storage!\n", mat->m, mat->n);
  for (ctr=0; ctr < nnz; ++ctr) {
    fscanf(fp, "%f ", &(mat->values.s[ctr]));
    if ( ctr % ((int)(nnz/100.0f))+1 == 0 ) { printf("val counter: %d\n",ctr); }
  }
  for (ctr=0; ctr < nnz; ++ctr) {
    fscanf(fp, "%d", &j);
    mat->rowind[ctr] = j;
    if ( ctr % ((int)(nnz/100.0f))+1 == 0 ) { printf("rowind counter: %d\n",ctr); }
  }
  for (ctr=0; ctr < mat->n+1; ++ctr) {
    fscanf(fp, "%d", &j);
    mat->colptr[ctr] = j;
    if ( ctr % ((int)(mat->n+1/100.0f))+1 == 0 ) { printf("colptr counter: %d\n",ctr); }
  }

  fclose(fp);
  return mat;
}

int
read_sparse( FILE *fp, double **vals, int *dim, int *cols )
{
  int i,j;
  double val;
  int fcode;
  int nnz;
  int checknnz;

  if (fscanf(fp, "SPARSE %d %d\n", dim, cols) != 2) { return -1; }
  if (fscanf(fp, "%d\n", &nnz) != 1) { return -1; }

  (*vals) = (double *)(calloc((*dim)*(*cols),sizeof(double)));

  if (*vals != NULL) { /* We can allocate this in one go. Do it the easy way. */

    for(checknnz=0;(fcode = fscanf(fp,"%d %d %lg",&i,&j,&val)) == 3;checknnz++) {

      i--; j--; /* Adjust for 1-based indexing */
      if (i < 0 || i > *dim-1 || j < 0 || j > *cols-1) { return -1; }
      (*vals)[i*(*cols) + j] = val;

    }

    /* We got here because the fscanf didn't work. */

    if (fcode == EOF) {

      //printf("tsnnls_test: Read matrix with %d (claimed %d) nonzeros.\n",
      //     checknnz,nnz);
      return 0;

    } else { return -1; }

  } else {

    printf("tsnnls_test: Couldn't allocate construction buffer for %d x %d matrix.\n",
	   *dim,*cols);
    exit(-1);

  }

}

/* This function reads a dense matrix from Octave as described above */

int
read_mat( FILE* fp, double **vals, int *dim, int *cols)

/* Returns -1 on failure, 0 on success. */

{
  *dim = -1;
  *cols=-1;
  int r, c;
  int startc;
  char linein[4096];

  /* First, we scan for lines starting with # */

  while ((startc = fgetc(fp)) != EOF) {

    if (startc == '#') { /* We've started a comment line */

      fgets(linein,sizeof(linein)-1,fp);
      sscanf(linein," rows: %d",dim);
      sscanf(linein," columns: %d",cols);

    } else if (!isspace(startc)) { /* We've run over into actual data. */

      ungetc(startc,fp); break;

    }

  }

  if (*dim == -1 || *cols == -1) {

    if (VERBOSITY >= 5) {

      printf("tsnnls_test: File format error. An Octave text file must contain"
	     "             two lines in the form\n"
	     "\n"
	     "             # rows: <m> \n"
	     "             # columns: <n> \n"
	     "\n"
	     "             in an initial block of lines starting with '#'.\n"
	     "\n"
	     "             This file parsed rows %d and cols %d.\n",
	   *dim,*cols);

    }

    return -1;

  }

  /* We have parsed the topmatter, and read the size of the matrix. */

  (*vals) = (double*)malloc(sizeof(double)*(*dim)*(*cols));

  for( r = 0; r<(*dim); r++ ) {
    for( c=0; c< (*cols); c++ ) {

      if (fscanf( fp, "%lf ", &((*vals)[r*(*cols)+c]) ) != 1) { return -1; }

    }
  }

  return 0;
}

#ifndef __DBL_EPSILON__
#define __DBL_EPSILON__ 2.2204460492503131e-16
#endif

static void
tsnnls_test(taucs_ccs_matrix *A,taucs_double *realx,taucs_double *b)
{
  int xItr;
  double err = 0;
  double expectedError=0;

  double residual;
  double ttime;
  double *x;

  int    i;

  for (i=0;i<nsolvers;i++) {

    if (i>0) { fprintf(stderr,"        "); }

    start_clock();

    /* Passing a 0.0 here means we will ALWAYS perform the error-reducing final
     * step with LSQR.
     */

    if (solvers[i] == pjv) {

      fprintf(stderr,"pjv     ");
      x = t_snnls_pjv(A, b, &residual, 0.0, 0 );

    } else if (solvers[i] == fallback) {

      fprintf(stderr,"fallback");
      x = t_snnls_fallback(A, b, &residual, 0.0, 0 ) ;

    } else if (solvers[i] == tsnnls) {

      fprintf(stderr,"block3  ");
      x = t_snnls(A, b, &residual, 0.0, 0);

    } else if (solvers[i] == spiv) {

      fprintf(stderr,"spiv    ");
      x = t_snnls_spiv(A,b,&residual,0.0,0,nconstrained);

    } else {

      printf("tsnnls_test: Illegal solver in tsnnls_test.\n");
      exit(1);

    }

    ttime = end_clock();

    fprintf(stderr, " %3d x %3d   %6f ", A->m, A->n, ttime );

    if( x == NULL ) {

      fprintf(stderr, "(Solver failed)           FAIL\n");
      exit_status = 1;

    } else {

      // compute relative error using ||(x*-x)||/||x||

      double* diff = (double*)malloc(sizeof(double)*A->n);

      for( xItr=0; xItr<A->n; xItr++ )
	diff[xItr] = realx[xItr] - x[xItr];

      int incX = {1};
      //double xnorm;
      //err = cblas_dnrm2( A->n, diff, 1 );
      //err /= cblas_dnrm2( A->n, realx, 1 );

      err = DNRM2_F77(&(A->n),diff,&incX);

      /* Expected absolute error given from theory is approx. cond^2*eps */
      expectedError = taucs_rcond(A);
      expectedError = 1.0/expectedError;

      expectedError *= expectedError;
      expectedError *= __DBL_EPSILON__;

      fprintf(stderr, "%8e %8e ", err, expectedError );

      if( err > expectedError ) {
	fprintf( stderr, "FAIL.\n" );
	exit_status = 1;
      } else {
	fprintf( stderr, "pass.\n");
      }

      free(diff);
      free(x);

    }

  }

  taucs_ccs_free(A);
  free(realx);
  free(b);
}

static void
tlsqr_test(taucs_ccs_matrix *A, taucs_double *realx, taucs_double *b)
{
  int xItr;
  double* x;
  double ttime;

  int i;

  for(i=0;i<nsolvers;i++) {

    start_clock();

    if (solvers[i] == tlsqr) {

      fprintf(stderr,"tlsqr   ");
      x = t_lsqr(A, b );

    } else if (solvers[i] == SOLlsqr) {

      fprintf(stderr,"SOLlsqr ");
      x = lsqrwrapper(A, b);

    } else {

      printf("tsnnls_test: Illegal solver in tlsqr_test.\n");
      exit(1);

    }

    ttime = end_clock();

    fprintf(stderr,"%3d x %3d  %6f  ",A->m,A->n,ttime);

    if( x == NULL ) {

      fprintf(stderr," (Could not compute) FAIL\n");
      exit_status = 1;

      taucs_ccs_free(A);
      free(b);
      free(realx);

    } else {

      // compute relative error using ||(x*-x)||/||x||

      double* diff = (double*)malloc(sizeof(double)*A->n);
      double expectedError;

      for( xItr=0; xItr<A->n; xItr++ )
	diff[xItr] = realx[xItr] - x[xItr];

      double err = 0;
      int incX = {1};

      //err = cblas_dnrm2( A->n, diff, 1 );
      //err /= cblas_dnrm2( A->n, realx, 1 );

      err = DNRM2_F77(&(A->n),diff,&incX);
      err /= DNRM2_F77(&(A->n),realx,&incX);

      /* Expected relative error given from theory is approx. cond^2*eps */
      expectedError = ((double)1.0/taucs_rcond(A));
      expectedError *= expectedError;
      expectedError *= __DBL_EPSILON__;

      printf( "%8e %8e  ", err, expectedError );

      if( err > expectedError ) {

	fprintf( stderr, "  FAIL.\n");
	exit_status = 1;

      } else {

	fprintf( stderr, "  pass.\n");

      }

      free(diff);
      free(x);
      taucs_ccs_free(A);
      free(b);
      free(realx);

    }

  }

}

/*
 * This function provides an example of using lsqr for the sparse matrix structure
 * WITHOUT going through our solver.
 */
extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );

static double*
lsqrwrapper( taucs_ccs_matrix* Af, double* b )
{
  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  double       *xf_raw;

  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, Af->m, Af->n );

  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = Af->m;
  lsqr_in->num_cols = Af->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = 0;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1/(10*sqrt(DBL_EPSILON));
  lsqr_in->max_iter = 4*lsqr_in->num_cols;
  //lsqr_in->num_rows + lsqr_in->num_cols + 1000;
  lsqr_in->lsqr_fp_out = NULL;
  for( bItr=0; bItr<Af->m; bItr++ )
    {
      lsqr_in->rhs_vec->elements[bItr] = b[bItr];
    }
  /* Here we set the initial solution vector guess, which is
   * a simple 1-vector. You might want to adjust this value for fine-tuning
   * t_snnls() for your application
   */
  for( bItr=0; bItr<Af->n; bItr++ )
    {
      lsqr_in->sol_vec->elements[bItr] = 1;
    }

  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = sparse_lsqr_mult;

  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, Af );

  /* copy into xf_raw */
  xf_raw = (double*)malloc(sizeof(double)*Af->n);
  for( bItr=0; bItr<Af->n; bItr++ ) // not really bItr here, but hey
    xf_raw[bItr] = lsqr_out->sol_vec->elements[bItr];

  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );

  return xf_raw;
}

int QUIET = 0;

int main( int argc, char* argv[] )
{

  FILE *outfile;
  double *Avals = NULL, *bvals = NULL, *xvals = NULL;
  int Adim,Acols,bdim,bcols,xdim,xcols;
  taucs_ccs_matrix *A;
  int tverbosity = 0;

  const char (**aname) = NULL;
  const char (**bname) = NULL;
  const char (**xname) = NULL;

  int nfiles;
  int i;

#ifdef WITH_ARGTABLE2

  struct arg_file *arg_Afile = arg_filen("A",NULL,"<A file>",1,64000,"matrix for problem in .mat or .sparse format");
  struct arg_file *arg_bfile = arg_filen("b",NULL,"<b file>",1,64000,"rhs vector b in .mat format");
  struct arg_file *arg_xfile = arg_filen("x",NULL,"<x file>",0,64000,"solution vector x in .mat format");
  struct arg_lit  *arg_tsnnls = arg_lit0(NULL,"tsnnls","solve with tsnnls");
  struct arg_lit  *arg_pjv = arg_lit0(NULL,"pjv","solve with reference solver");
  struct arg_lit  *arg_tlsqr = arg_lit0(NULL,"tlsqr","solve with tlsqr");
  struct arg_lit  *arg_lsqr = arg_lit0(NULL,"lsqr","solve with SOL lsqr");
  struct arg_lit  *arg_fallback = arg_lit0(NULL,"fallback","solve with fallback");
  struct arg_lit  *arg_spiv = arg_lit0(NULL,"spiv","solve with single pivoting solver");
  struct arg_int  *arg_constrained = arg_int0(NULL,"nc","<0-m>","number of constrained vars (spiv only)");

  struct arg_int  *arg_verb = arg_int0("v","Verbosity","<0-10>","verbosity for tsnnls solver");

  struct arg_lit  *arg_help = arg_lit0("?","help","display help message");

  struct arg_end *end = arg_end(20);

  void *argtable[] = {arg_Afile,arg_bfile,arg_xfile,arg_tsnnls,
		      arg_pjv,arg_fallback,arg_tlsqr,arg_lsqr,arg_spiv,
		      arg_constrained,arg_verb,arg_help,end};

  int nerrors;

  /* Now we parse the arguments */

  if (arg_nullcheck(argtable) != 0) {
    /* NULL entries detected, allocations must have failed. */
    fprintf(stderr,"%s: insufficient memory\n",argv[0]);
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  nerrors = arg_parse(argc,argv,argtable);

  /* special case: '--help' takes precedence over error reporting */
  if (arg_help->count > 0) {

    printf("tsnnls_test solves a least-squares or constrained least-squares\n"
	   "problem in the A x = b. The files are expected to be in a sparse\n"
	   "format or in the Octave text format.\n"
	   "\n");

    fprintf(stderr,"Usage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

    printf("If a solution x is provided, tsnnls_test compares the results of the\n"
	   "calculation with the given solution, provides timing information and\n"
	   "returns 0 if the results match the given solution and 1 otherwise.\n"
	   "\n"
	   "If no solution is provided, tsnnls_test writes the solution to x.mat.\n");

    return 0;
  }

  /* If the parser returned any errors then display them and exit */
  if (nerrors > 0) {
    /* Display the error details contained in the arg_end struct.*/
    fprintf(stderr,"\n");
    arg_print_errors(stderr,end,argv[0]);
    fprintf(stderr,"\nUsage: %s ",argv[0]);
    arg_print_syntax(stderr,argtable,"\n");
    arg_print_glossary(stderr,argtable,"  %-25s %s\n");
    arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));
    return 1;
  }

  /* Now we read the arguments and move into local vars. */

  if (arg_Afile->count != arg_bfile->count) {

    printf("tsnnls_test: Input %d A matrices, but %d b right-hand vectors.\n"
	   "             These numbers must be the same to run multiple files.\n",
	   arg_Afile->count,arg_bfile->count);
    exit(1);

  }

  nfiles = arg_Afile->count;

  aname = calloc(sizeof(char *),nfiles);
  bname = calloc(sizeof(char *),nfiles);

  for(i=0;i<nfiles;i++) {

    aname[i] = arg_Afile->filename[i];
    bname[i] = arg_bfile->filename[i];

  }

  if (arg_xfile->count > 0) {

    if (arg_xfile->count != nfiles) {

      printf("tsnnls: Number of solutions x (%d) must match number of A and b files (%d).\n",
	     arg_xfile->count,nfiles);

      exit(1);

    }

    xname = calloc(sizeof(char *),nfiles);
    for(i=0;i<nfiles;i++) { xname[i] = arg_xfile->filename[i]; }

  }

  if (arg_lsqr->count > 0)  { solvers[nsolvers++] = SOLlsqr; }
  if (arg_tsnnls->count > 0) { solvers[nsolvers++] = tsnnls; }
  if (arg_pjv->count > 0)    { solvers[nsolvers++] = pjv; }
  if (arg_tlsqr->count > 0)  { solvers[nsolvers++]  = tlsqr; }
  if (arg_fallback->count > 0) {solvers[nsolvers++] = fallback; }
  if (arg_spiv->count > 0) { solvers[nsolvers++] = spiv; }
  if (arg_verb->count > 0) { tverbosity = arg_verb->ival[0]; }

  if (arg_constrained->count > 0) {nconstrained = arg_constrained->ival[0];}

  arg_freetable(argtable,sizeof(argtable)/sizeof(argtable[0]));

#else

  if (argc < 4 || (strcmp(argv[argc-1],"--tsnnls") && strcmp(argv[argc-1],"--tlsqr") && \
		   strcmp(argv[argc-1],"--spiv") && strcmp(argv[argc-1],"--pjv") && \
		   strcmp(argv[argc-1],"--fallback"))
      || argc > 6) {

    printf("Usage: tsnnls_test <A file> <b file> "
	   "<x file (optional)> <--fallback|--tsnnls|--tlsqr|--pjv|--spiv>\n");
    exit(1);

  }

  if (!strcmp(argv[1],"--help") || !strcmp(argv[1],"-?")) {

    printf("tsnnls_test solves a least-squares or constrained least-squares\n"
	   "problem in the A x = b. The files are expected to be in a sparse\n"
	   "format or in the Octave text format.\n"
	   "\n");

    printf("Usage: tsnnls_test <A file> <b file> <x file (optional)> "
	   "<--tsnnls|--tlsqr|--pjv|--spiv>\n\n");

    printf("If a solution x is provided, tsnnls_test compares the results of the\n"
	   "calculation with the given solution, provides timing information and\n"
	   "returns 0 if the results match the given solution and 1 otherwise.\n"
	   "\n"
	   "If no solution is provided, tsnnls_test writes the solution to x.mat.\n");

    exit(1);

  }

  /* Now we attempt to parse the arguments and load files. */

  nfiles = 1;

  aname = (const char**)calloc(sizeof(char *),1);
  bname = (const char**)calloc(sizeof(char *),1);

  aname[0] = argv[1];
  bname[0] = argv[2];

  if (argc == 5) {

    xname = (const char**)calloc(sizeof(char *),1);
    xname[0] = argv[3];

  }

  nsolvers=1;

  if (!strcmp(argv[argc-1],"--tsnnls")) { solvers[0] = tsnnls; }
  else if (!strcmp(argv[argc-1],"--pjv")) { solvers[0] = pjv; }
  else if (!strcmp(argv[argc-1],"--tlsqr")) { solvers[0] = tlsqr; }
  else if (!strcmp(argv[argc-1],"--fallback")) { solvers[0] = fallback; }
  else if (!strcmp(argv[argc-1],"--lsqr")) { solvers[0] = SOLlsqr; }
  else if (!strcmp(argv[argc-1],"--spiv")) { solvers[0] = spiv; }
  else {

    printf("tsnnls_test: Unknown solver %s.\n",argv[argc-1]);
    exit(1);

  }

#endif

  printf("tsnnls_test %s (%s %s).\n",PACKAGE_VERSION, __DATE__ , __TIME__ );
  tsnnls_version(NULL,0);
  tsnnls_verbosity(tverbosity);

  fprintf(stderr,"Loaded  Solver   matrixdims  runtime  error        errbound     Result\n");
  fprintf(stderr,"----------------------------------------------------------------------\n");

  int flags = TAUCS_SINGLE;
  for(i=0;i<nfiles;i++) {
    printf("Reading %s\n",aname[i]);
    A = read_sparse_non_stupid( aname[i], flags );
    Adim = A->m;
    Acols = A->n;
    //Avals = loadvals(aname[i], &Adim, &Acols);
    //A = taucs_construct_sorted_ccs_matrix(Avals, Acols, Adim);
    //free(Avals);

    fprintf(stderr,"A");

    if (nconstrained == -1) { nconstrained = Acols; }

    /* Now load b. */

    bvals = loadvals(bname[i], &bdim, &bcols);

    if (bdim != Adim || bcols != 1) {

      printf("tsnnls_test: We expect the second (b) file %s to be %d x 1.\n"
	     "             The given file is %d x %d.\n",
	     bname[i],Adim,bdim,bcols);
      exit(1);

    }

    fprintf(stderr,"b");

    /* Now, if present, load x */

    if (xname != NULL) {

      xvals = loadvals(xname[i], &xdim, &xcols);

      if (xdim != Acols || xcols != 1) {

	printf("tsnnls_test: We expect the third (x) file %s to be %d x 1.\n"
	       "             The given file is %d x %d.\n",
	       xname[i],Acols,xdim,xcols);
	exit(1);

      }

      fprintf(stderr,"x     ");

      int solItr;

      if (solvers[0] == tsnnls || solvers[0] == pjv ||
	  solvers[0] == fallback || solvers[0] == spiv) {

	for(solItr = 0;solItr < nsolvers;solItr++) {

	  if (solvers[solItr] == tlsqr || solvers[solItr] == SOLlsqr) {

	    fprintf(stderr,
		    "\ntsnnls: Error. Can't mix least squares solvers (tlsqr, SOLlsqr)\n"
		    "         with constrained solvers (tsnnls, pjv, spiv, fallback).\n");
	    exit(1);

	  }

	}

	tsnnls_test(A,xvals,bvals);

      } else if (solvers[0] == tlsqr || solvers[0] == SOLlsqr) {

	for (solItr = 0;solItr < nsolvers;solItr++) {

	  if ((solvers[solItr] != tlsqr) && (solvers[solItr] != SOLlsqr)) {

	    fprintf(stderr,
		    "\ntsnnls: Error. Can't mix least squares solvers (tlsqr, SOLlsqr)\n"
		    "         with constrained solvers (tsnnls, pjv, spiv, fallback).\n");
	    exit(1);

	  }

	}

	tlsqr_test(A,xvals,bvals);

      }

    } else {

      /* If we've survived this long, we're in problem mode instead of test mode. */
      fprintf(stderr,"      ");
      double residual;

      start_clock();

      if (solvers[0] == fallback) {

	fprintf(stderr,"fallback");
	xvals = t_snnls_fallback(A, bvals, &residual, -1, 1);

      } else if (solvers[0] == pjv) {

	fprintf(stderr,"pjv     ");
	xvals = t_snnls_pjv(A, bvals, &residual, -1, 1);

      } else if (solvers[0] == spiv) {

	fprintf(stderr,"spiv    ");
	xvals = t_snnls_spiv(A, bvals, &residual, -1, 1,nconstrained);

      } else if (solvers[0] == tlsqr) {

	fprintf(stderr,"tlsqr   ");
	xvals = t_lsqr(A, bvals);

      } else if (solvers[0] == SOLlsqr) {

	fprintf(stderr,"SOLlsqr ");
	xvals = lsqrwrapper(A, bvals);

      } else {

	fprintf(stderr,"block3  ");
	xvals = t_snnls(A, bvals, &residual, -1, 1);

      }

      double ttime;
      ttime = end_clock();

      fprintf(stderr," %3d x %3d",A->m,A->n);
      fprintf(stderr,"  %6f ",ttime);

      if (xvals == NULL) {

	fprintf(stderr," Fail.\n");
	exit_status = 1;

      } else {

	printf(" Success.\n");

	char outfilename[256];

	sprintf(outfilename,"x%d.mat",i);
	outfile = fopen(outfilename,"w");

	if (outfile == NULL) {

	  printf("tsnnls_test: Could not open %s to write solution.\n",outfilename);
	  exit(1);

	}

	colvector_write_mat(outfile,xvals,Acols,"x");
	free(xvals);

      }

      free(bvals);
      fclose(outfile);
      taucs_ccs_free(A);

    }

  }

  free(aname);
  free(bname);
  if (xname != NULL) { free(xname); }

  exit(exit_status);

}

int old_readsparse(FILE *fp,double **x,double **b,taucs_ccs_matrix **outMatrix)

{
  int dim = 0, cols=0, r, nnz, colOffset, cItr;
  double* vals = NULL;
  taucs_ccs_matrix* result;
  int vItr=0;
  int*	colcounts;
  double** colels;
  int**	colrows;
  int*	colitrs;
  fpos_t  elementsStart;
  int elr, elc;
  double theVal;

  fgetpos(fp, &elementsStart);

  if (fscanf(fp, "SPARSE %d %d\n", &dim, &cols) != 2) { return -1; }
  if (fscanf(fp, "%d\n", &nnz) != 1) { return -1; }

  result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  result->n = cols;
  result->m = dim;
  result->flags = TAUCS_DOUBLE;

  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  result->rowind = (int*)malloc(sizeof(int)*nnz);
  result->values.d = (double*)malloc(sizeof(double)*nnz);

  vals = (double*)calloc(nnz,sizeof(double));
  colels = (double**)calloc(result->n,sizeof(double*));
  colrows = (int**)calloc(result->n,sizeof(int*));
  colcounts = (int*)calloc(result->n,sizeof(int));

  /* we have to munger our file reading a bit since we're given
   * output in row1 col1, col2, ... row2 col1, col2, but it's better
   * than allocating a m*n double array when m and n are large
   */
  for( r=0; r<nnz; r++ )
    {
      /* first get column counts */
      fscanf(fp, "%d %d %lf\n", &elr, &elc, &theVal);
      elr--;
      elc--; // adjust for 0 based indexing
      colcounts[elc]++;
    }

  /* Now that we have the count of each column, we can allocate colels to
   * store the values on our second pass
   */
  for( cItr=0; cItr<result->n; cItr++ )
    {
      colels[cItr] = (double*)calloc(colcounts[cItr], sizeof(double));
      colrows[cItr] = (int*)calloc(colcounts[cItr], sizeof(int));
    }

  fsetpos(fp, &elementsStart);

  colitrs = (int*)calloc(result->n, sizeof(int));
  for( r=0; r<nnz; r++ )
    {
      fscanf(fp, "%d %d %lf\n", &elr, &elc, &theVal);
      elr--;
      elc--; // adjust for 0 based indexing
      colels[elc][colitrs[elc]] = theVal;
      colrows[elc][colitrs[elc]] = elr;
      colitrs[elc]++;
    }
  free(colitrs);

  /* now that we have the elements ordered in their columns and we've also
   * kept track of the row indices, we are ready to construct the ccs structure
   * which was (comparatively easy!) to allocate
   */
  colOffset = 0;
  for( cItr=0; cItr<result->n; cItr++ )
    {
      result->colptr[cItr] = colOffset;
      for( vItr=0; vItr<colcounts[cItr]; vItr++ )
	{
	  result->values.d[colOffset] = colels[cItr][vItr];
	  result->rowind[colOffset] = colrows[cItr][vItr];
	  colOffset++;
	}
    }
  result->colptr[result->n] = nnz;

  free(vals);
  free(colcounts);
  for( cItr=0; cItr<result->n; cItr++ )
    free(colels[cItr]);
  for( cItr=0; cItr<result->n; cItr++ )
    free(colrows[cItr]);
  free(colels);
  free(colrows);

  *x = (double*)malloc(sizeof(double)*cols);
  for( r=0; r<cols; r++ )
    {
      double a;
      fscanf(fp, "%lf\n", &a);
      (*x)[r] = a;
    }

  *b = (double*)malloc(sizeof(double)*dim);
  for( r=0; r<dim; r++ )
    fscanf(fp, "%lf\n", &(*b)[r]);

  *outMatrix = result;

  return 0;
}

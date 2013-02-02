
/*
 * This program is free software distributed under the LGPL. A copy of
 * the license should have been included with this archive in a file
 * named 'LICENSE'. You can read the license there or on the web at:
 * http://www.gnu.org/licenses/lgpl-3.0.txt
 */

#include <config.h>
#include "tsnnls_blas_wrappers.h"

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#ifdef HAVE_FLOAT_H
  #include <float.h>
#endif 

//#ifdef HAVE_DARWIN  /* We expect to use the Accelerate framework */
//  #include <vecLib/vBLAS.h>
//  #include <vecLib/clapack.h>
//#else
//  #include "gsl_cblas.h"
//#endif 

#ifdef HAVE_CLAPACK_H
  #include <clapack.h>
#else 
  #ifdef HAVE_ATLAS_CLAPACK_H
     #include <atlas/clapack.h>
  #else
     #ifdef HAVE_VECLIB_CLAPACK_H
       #include <vecLib/clapack.h>
     #endif
  #endif
#endif

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STRING_H
  #include <string.h>
#endif

#include "acint32_type.h"
#include "lsqr.h"
#include "tsnnls.h"

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

// 

/* This file contains code for the TAUCS version of lsqr, including
   two internal procedures called by the main routine below. 

   Cantarella/Piatek, 5/04. 
*/

/* Procedure is internal to taucs_ccs_aprime_times_a, which is 
 * presented below. 
 */
double 
taucs_dotcols( const taucs_ccs_matrix* A, int col1, int col2 )
{
	double val = 0;
	int i, j;

	i=A->colptr[col1];
	j=A->colptr[col2];

	// eliminate some loop invariant pointer derefs
	double* Avals = A->values.d;
	int* rowinds = A->rowind;
	int* colptrs = A->colptr;
	
	int stopI, stopJ;
	
	stopI = colptrs[col1+1];
	stopJ = colptrs[col2+1];

	while( 1==1 )
	{
		if( rowinds[i] == rowinds[j] )
			val += Avals[i++]*Avals[j++];
		else if( rowinds[i] < rowinds[j] )
			i++;
		else
			j++;
		if( i>=stopI )
			break;
		if( j>=stopJ )
			break;
	}
	return val;
}

/* 
 * This routine computes A'*A for full A. This is useful for comparing
 * sparse/nonsparse performance to determine if your problem can
 * benefit from a different representation.  It is left in the archive
 * as a convenience
 */
double* 
full_aprime_times_a( double* A, int rows, int cols )
{
  int rItr, cItr, colOffset;
  double* result = (double*)calloc(cols*cols, sizeof(double));
  int incX,incY,N;

  colOffset = 0;
  N = rows; 
  incX = incY = cols;
  
  for( cItr=0; cItr<cols; cItr++ )  {
    for( rItr=cItr; rItr<cols; rItr++ ) {
      //result[rItr*cols + cItr] = cblas_ddot( rows, &A[cItr], 
      //                                       cols, &A[rItr], cols );
      result[rItr*cols + cItr] = DDOT_F77(&N, &A[cItr], &incX, 
			                  &A[rItr], &incY);
    }
  }
  
  return result;
}

/* 
 * Procedure computes the matrix product A'*A for a sparse matrix 
 * A in compressed column storage. The result is a symmetric TAUCS
 * matrix, which is also in CCS format. 
 */
taucs_ccs_matrix* 
taucs_ccs_aprime_times_a( taucs_ccs_matrix* A )
{
  int rItr, cItr, colOffset; // maxSize, scItr, srItr;
  taucs_ccs_matrix* result = (taucs_ccs_matrix*)malloc(sizeof(taucs_ccs_matrix));
  register double v;
  int newSize, currentSize;

  result->m = A->n;
  result->n = A->n;
  
  /* we'll be symmetric */
  result->flags = TAUCS_DOUBLE;
  result->flags = result->flags | TAUCS_SYMMETRIC;
  result->flags = result->flags | TAUCS_LOWER; // rep the lower half
  
  result->colptr = (int*)malloc(sizeof(int)*(result->n+1));
  
  currentSize = A->colptr[A->n]*2; // start with the number of entries in A*2 and see if we need more
  result->values.d = (double*)calloc(sizeof(taucs_double),currentSize); // Changed to track down UIval bug
  result->rowind = (int*)malloc(sizeof(int)*currentSize);
  
  colOffset = 0;
  double* valsPtr = result->values.d;
  int* rowptrs = result->rowind;
  int* colptrs = result->colptr;
  int Acols = A->n;
  //  int* Acolptrs = A->colptr;
  // double* AvalsPtr = A->values.d;
  // int* Arowptrs = A->rowind;
  // int stopPoint, interiorStopPoint;
  
  for( cItr=0; cItr<Acols; cItr++ ) {
    
    colptrs[cItr] = colOffset;
    
    for( rItr=cItr; rItr<Acols; rItr++ ) {
      v = taucs_dotcols(A,cItr,rItr);
      
      if( v == 0.0 ) {
	continue;
      }
      else {
	valsPtr[colOffset] = v;
	rowptrs[colOffset] = rItr;
	colOffset++;
	
	if( colOffset < currentSize )
	  continue;
	else {
	  /* we need to increase our allocation size.  */
	  newSize = 2*currentSize;
	  int* newRows = (int*)realloc(rowptrs, sizeof(int)*newSize);
	  double* newVals = (double*)realloc(valsPtr, sizeof(double)*newSize);
	  
	  currentSize = newSize;
	  
	  if( newRows == NULL || newVals == NULL ) {
	    fprintf( stderr, "tsnnls: Out of memory!\n" );
	  }
	  
	  result->values.d = newVals;
	  valsPtr = newVals;
	  result->rowind = newRows;
	  rowptrs = newRows;
	}
      }
    }
  }
  colptrs[cItr] = colOffset; 
  return result;
}

void
ccs_to_lapack( taucs_ccs_matrix* L, double** lapackL, int* N, int* LDA, double* ANORM )
{	
  /* Construct LAPACK representation of A and compute the 1 norm of A */
  int vSize;
  int cItr, rItr;
  int rowCount = L->m;
  double localMax = 0;
  
  *ANORM = 0;
  
  if( (L->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
    {
      vSize = L->n*L->n;
      rowCount = L->n;
    }
  else
    vSize = L->m*L->n;
	
  *lapackL = (double*)calloc(vSize,sizeof(double));
  
  for( cItr=0; cItr<L->n; cItr++ )
    {
      localMax = 0;
      for( rItr=L->colptr[cItr]; rItr<L->colptr[cItr+1]; rItr++ )
	{
	  int index = -1;
	  index = L->rowind[rItr] + cItr*rowCount;
	  (*lapackL)[index] = L->values.d[rItr];
	  localMax += fabs(L->values.d[rItr]);
	}
      if( localMax > *ANORM )
	*ANORM = localMax;
    }
  
  *N = L->n;
  *LDA = rowCount;
}

static double
t_condest( void* mfR )
{
  #ifndef HAVE_ATLAS_LAPACK 

  /* If we have a full LAPACK available, use dpocon to estimate condition number. */
  
  taucs_ccs_matrix* L;
  double* lapackL;
  ACINT32_TYPE N, LDA, INFO;  /* Lapack expects 32 bit ints as integers. */
  char	UPLO;
  double  ANORM = 0;
  double* WORK;
  ACINT32_TYPE *IWORK;
  double  RCOND;
  
  L = taucs_supernodal_factor_to_ccs(mfR);
  
  /* Construct LAPACK representation of A and compute the 1 norm of A */
  int vSize;
  int cItr, rItr;
  int rowCount = L->m;
  double localMax = 0;
  
  if( (L->flags & TAUCS_SYMMETRIC)==TAUCS_SYMMETRIC )
    {
      vSize = L->n*L->n;
      rowCount = L->n;
	}
  else
    vSize = L->m*L->n;
  
  lapackL = (double*)calloc(vSize,sizeof(double));
  
  for( cItr=0; cItr<L->n; cItr++ )
    {
      localMax = 0;
      for( rItr=L->colptr[cItr]; rItr<L->colptr[cItr+1]; rItr++ )
	{
	  int index = -1;
	  index = L->rowind[rItr] + cItr*rowCount;
	  lapackL[index] = L->values.d[rItr];
	  localMax += fabs(L->values.d[rItr]);
	}
      if( localMax > ANORM )
	ANORM = localMax;
    }
  
  N = L->n;
  LDA = L->m;
  UPLO = 'L';
  WORK = (double*)malloc(sizeof(double)*3*N);
  IWORK = (ACINT32_TYPE*)malloc(sizeof(ACINT32_TYPE)*N);
  
  dpocon_( &UPLO, &N, lapackL, &LDA, &ANORM, &RCOND, WORK, IWORK, &INFO );
  
  free(WORK);
  free(IWORK);
  taucs_ccs_free(L);
  free(lapackL);
  
  return RCOND;

  #else 

  /* We have only a limited ATLAS LAPACK available, which doesn't include the condition
     number estimating code. We can still get a condition number estimate from lsqr. */
 
  taucs_ccs_matrix *A;

  lsqr_input   *lsqr_in;
  lsqr_output  *lsqr_out;
  lsqr_work    *lsqr_work;
  lsqr_func    *lsqr_func;
  int bItr;
  
  double        cond;

  extern void sparse_lsqr_mult( long mode, dvec* x, dvec* y, void* prod );

  A = taucs_supernodal_factor_to_ccs(mfR);  
  alloc_lsqr_mem( &lsqr_in, &lsqr_out, &lsqr_work, &lsqr_func, A->m, A->n );
  
  /* we let lsqr() itself handle the 0 values in this structure */
  lsqr_in->num_rows = A->m;
  lsqr_in->num_cols = A->n;
  lsqr_in->damp_val = 0;
  lsqr_in->rel_mat_err = 0;
  lsqr_in->rel_rhs_err = 0;
  lsqr_in->cond_lim = 1e16;
  lsqr_in->max_iter = lsqr_in->num_rows + lsqr_in->num_cols + 1000;
  lsqr_in->lsqr_fp_out = NULL;	
  
  for( bItr=0; bItr<A->m; bItr++ )
    {
      lsqr_in->rhs_vec->elements[bItr] = 1;  /* We will solve Ax = [1 1 .... 1]^T. */
    }
  /* Here we set the initial solution vector guess, which is 
   * a simple 1-vector. You might want to adjust this value for fine-tuning
   * t_snnls() for your application
   */
  for( bItr=0; bItr<A->n; bItr++ ) {
    lsqr_in->sol_vec->elements[bItr] = 1; 
  }
  
  /* This is a function pointer to the matrix-vector multiplier */
  lsqr_func->mat_vec_prod = sparse_lsqr_mult;
  
  lsqr( lsqr_in, lsqr_out, lsqr_work, lsqr_func, A );  
  cond = lsqr_out->mat_cond_num;
  
  free_lsqr_mem( lsqr_in, lsqr_out, lsqr_work, lsqr_func );
  taucs_ccs_free(A);
  
  return cond;

  #endif
}

/*
 * This is a custom version of the TAUCS based lsqr solver for use
 * with the t_snnls() block principle pivoting function, which
 * precomputes several values of interest to lsqr so as to not repeat
 * their calcuation during the multiple pivots which can occur.
 */
taucs_double*
t_snnlslsqr(taucs_ccs_matrix *A,
	    taucs_double *b, 
	    taucs_ccs_matrix* ApA, 
	    int* F,
	    double* outRcond )
{
  taucs_ccs_matrix   *ApAperm;
  taucs_double *Apb /*n x 1*/, \
    *ApAx /*n x 1*/, \
    *x /* n x 1 */, \
    *Itstep /*n x 1*/; 
  double	   *x_unscrambled;
  void		   *mfR; /* The Cholesky factor R, stored by TAUCS. */
  int		   *perm, *invperm;
  char		   ordering[1024] = "amd";
  
  //ordering = getenv("COL_ORDERING");
  //if( ordering == NULL )
  //   {
      /* use amd ordering if the user hasn't specified anything else */
  //    putenv("COL_ORDERING=amd");
  //    ordering = getenv("COL_ORDERING");
  //  }
  taucs_ccs_order(ApA, &perm, &invperm, ordering);
  ApAperm = taucs_ccs_permute_symmetrically(ApA, perm, invperm);
  
  Apb = (taucs_double*)calloc(A->m, sizeof(taucs_double));
  
  mfR = taucs_ccs_factor_llt_mf(ApAperm);
  if( mfR == NULL )
    {
      // free the memory we've allocated so far and return
      taucs_ccs_free(ApAperm);
      free(Apb);
      free(perm);
      free(invperm);
      return NULL;
    }
  
  if( outRcond != NULL )
    *outRcond = t_condest(mfR);
  
  /* We now solve the first equation: x = R\(A'*A*b). */
  x = (taucs_double*)calloc(A->n,sizeof(taucs_double));
  
  taucs_transpose_vec_times_matrix_nosub(b, A, Apb);
  // we have to permute A'b to be meaningful.
  {
    double* apbperm = Apb;
    Apb = (double*)malloc(sizeof(double)*A->n);
    taucs_vec_permute(A->n, TAUCS_DOUBLE, apbperm, Apb, perm);
    free(apbperm);
  }
  
  taucs_supernodal_solve_llt(mfR,x,Apb);  /* n x n * n x 1 = n x 1 */
  
  /* Given the base solution, we now update it by refinement. */
  ApAx = (taucs_double*)malloc(sizeof(double)*A->n);
  Itstep = (taucs_double*)calloc(A->n,sizeof(taucs_double)); 
	
  double* scratch = (double*)malloc(sizeof(double)*ApAperm->n);
  memcpy(scratch, Apb, sizeof(double)*ApAperm->n);
  
  ourtaucs_ccs_times_vec(ApAperm,x,ApAx);  /* n x n * n x 1 = n * 1. */

  double alpha = {-1.0};
  int    incX = {1},incY = {1};

  //cblas_daxpy(A->n,-1.0,(double *)(ApAx),1,(double *)(scratch),1); 
  /* Apb = Apb - ApAx */
  DAXPY_F77(&(A->n),&alpha,(double *)(ApAx),&incX,(double *)(scratch),&incY);

  //refinementEps = cblas_dnrm2(A->n, scratch, 1);
  taucs_supernodal_solve_llt(mfR,Itstep,scratch); /* Itstep = R\Apb. */

  alpha = 1.0;
  //cblas_daxpy(A->n,1.0,(double *)(Itstep),1,(double *)(x),1);  
  DAXPY_F77(&(A->n),&alpha,(double *)(Itstep),&incX,(double *)(x),&incY);
  /* x = x + Itstep */
  
  free(scratch);
  free(Itstep);
  free(ApAx);
  
  /* We now have a solution, but must clean up the memory we allocated
   * and unscramble x 
   */
  x_unscrambled = (double*)calloc(sizeof(double),ApA->n);
  taucs_vec_permute(ApA->n, TAUCS_DOUBLE, x, x_unscrambled, invperm);
  
  taucs_ccs_free(ApAperm);
  free(Apb); 
  free(perm);
  free(invperm);
  free(x);
  taucs_supernodal_factor_free(mfR);
  
  return x_unscrambled;
}

taucs_double*
t_lsqr(taucs_ccs_matrix *A, taucs_double *b)
{
  int *F;
  int i;
  double* x;
  taucs_ccs_matrix* apda;

  F = (int*)malloc(sizeof(int)*A->n);
  for( i=0; i<A->n; i++ )
    F[i] = i;
  apda = taucs_ccs_aprime_times_a(A);
  x = t_snnlslsqr(A, b, apda, F, NULL);
  taucs_ccs_free(apda);
  free(F);
  
  return x;
}

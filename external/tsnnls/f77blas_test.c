/* 

f77blas_test.c 

An f77blas test program. We generate a random matrix, copy it, and
perform a sequence of cblas operations which should take us back to
the original matrix. The objective is to make sure that the cblas
installation on this system is functional.

*/

#include <config.h>
#include "tsnnls_blas_wrappers.h"

#ifdef HAVE_MATH_H
  #include <math.h>
#endif

#ifdef HAVE_STDLIB_H
  #include <stdlib.h>
#endif

#ifdef HAVE_STDIO_H
  #include <stdio.h>
#endif

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#ifdef __cplusplus
extern "C" /* Avoid C++ name mangling */
#endif

//#define DAXPY_F77 F77_FUNC (daxpy,DAXPY)
//#define DSYMV_F77 F77_FUNC (dsymv,DSYMV)
//#define DSYRK_F77 F77_FUNC (dsyrk,DSYRK)

#define VECLEN 1024
#define MATDIM 73

int main()
{

  double vecA[VECLEN], vecB[VECLEN], vecC[VECLEN];
  int i,incA={1},incB={1},N={VECLEN};
  double alpha = -2.3;

  fprintf(stderr,"f77blas tests\n");

  /* L1 BLAS test for "daxpy" */

  fprintf(stderr,"%s test... ","DAXPY_F77");

  for(i=0;i<VECLEN;i++) {

    vecC[i] = vecA[i] = drand48();
    vecB[i] = drand48();

  }

  DAXPY_F77(&N,&alpha,vecB,&incB,vecA,&incA);  /* vecA += alpha*vecB */

  for(i=0;i<VECLEN;i++) {
    
    if (fabs(vecC[i] + alpha*vecB[i] - vecA[i]) > 1e-12) {

      fprintf(stderr," failed.\n"
	      "A[%d] should have been %g * %g + %g = %g, was %g.\n\n",
	      i,alpha,vecB[i],vecC[i],alpha*vecB[i]+vecC[i],vecA[i]);

      exit(1);
    }

  }

  fprintf(stderr," passed.\n");

  /* L2 BLAS check for dsymv */

  double vecX[MATDIM],vecY[MATDIM],A[MATDIM*MATDIM],\
    vecYsto[MATDIM],vecCheck[MATDIM];
  int j,lda = {MATDIM},incX={1},incY={1};
  char uplo = {'U'};
  double beta = 1.66;

  fprintf(stderr,"%s test... ","DSYMV_F77");

  N = MATDIM; /* Reset N to dimension of NxN matrix A */

  for(i=0;i<MATDIM;i++) {

    vecX[i] = drand48();
    vecY[i] = vecYsto[i] = drand48();

    for(j=i;j<MATDIM;j++) {

      /* A(j,i) = A(i,j) = drand48();
	 We used FORTRAN column-major order for the matrix here. */
      
      A[j + i*MATDIM] = A[i + j*MATDIM] = drand48();  
      
    }

    vecCheck[i] = 0;

  }

  DSYMV_F77(&uplo,&N,&alpha,A,&lda,vecX,&incX,&beta,vecY,&incY);

  /* Should compute vecY = alpha A vecX + beta vecY */

  for(i=0;i<MATDIM;i++) {

    for(j=0;j<MATDIM;j++) {

      vecCheck[i] += alpha*(A[i + j*MATDIM]*vecX[j]);

    }

    vecCheck[i] += beta*vecYsto[i];

  }

  for(i=0;i<MATDIM;i++) {

    if (fabs(vecCheck[i] - vecY[i]) > 1e-12) {

      fprintf(stderr," failed.\n"
	      "vecCheck[%d] was %g and vecY[%d] was %g (should be equal).\n\n",
	      i,vecCheck[i],i,vecY[i]);
      exit(1);

    }

  }

  fprintf(stderr," passed.\n");

  /* L3 BLAS test. dsyrk. */
  /* This tests an NxK matrix newA, where N = MATDIM, K = 2*MATDIM. */

  double newA[MATDIM*2*MATDIM],C[MATDIM*MATDIM], \
    newAnewAT[MATDIM*MATDIM],Csto[MATDIM*MATDIM];
  char trans = 'N';
  int k;
  int K = {2*MATDIM};

  fprintf(stderr,"%s test... ","DSYRK_F77");

  for(i=0;i<MATDIM;i++) {

    for(j=0;j<MATDIM;j++) {

      newAnewAT[i + MATDIM*j] = 0;

    }

    for(j=i;j<MATDIM;j++) { /* C is a MATDIM*MATDIM symmetric matrix */

      Csto[i + MATDIM*j] = Csto[j + MATDIM*i] = 
	C[i + MATDIM*j] = C[j + MATDIM*i] = drand48();

    }

    for(j=0;j<2*MATDIM;j++) {

      newA[i + MATDIM*j] = drand48();

    }

  }

  DSYRK_F77(&uplo,&trans,&N,&K,&alpha,newA,&N,&beta,C,&N);

  /* Computes C = alpha newAnewA^T + beta C, which is a MATDIM*MATDIM matrix. */

  for(i=0;i<MATDIM;i++) {

    for(j=0;j<MATDIM;j++) {

      for(k=0;k<2*MATDIM;k++) {

	newAnewAT[i + MATDIM*j] += newA[i + MATDIM*k] * newA[j + MATDIM*k]; /* newAnewA^T(i,j) */

      }

      newAnewAT[i + MATDIM*j] *= alpha;
      newAnewAT[i + MATDIM*j] += beta*Csto[i + MATDIM*j];

    }

  }

  for(i=0;i<MATDIM;i++) {

    for(j=i;j<MATDIM;j++) {  /* Only the upper triangle is set by dsyrk */

      if (fabs(newAnewAT[i + MATDIM*j] - C[i + MATDIM*j]) > 1e-12) {

	fprintf(stderr," failed.\n"
		"C(%d,%d) was %g, should have been %g.\n",
		i,j,C[i+MATDIM*j],newAnewAT[i+MATDIM*j]);

	exit(1);

      }

    }

  }
  
  fprintf(stderr," passed.\n");

  /* End of tests */

  exit(0);

}
  

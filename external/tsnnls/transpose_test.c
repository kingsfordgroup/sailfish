/*

  transpose_test.c: This file runs a quick test of the
  taucs_ccs_transpose function by generating a bunch of random
  matrices and taking their transpose to make sure that the system
  works.

*/

#define NUM_TESTS 20
#define Msize 40
#define Nsize 555

#include<config.h>

#include"acint32_type.h"
#include"tsnnls_blas_wrappers.h"

#ifdef HAVE_STDLIB_H
 #include<stdlib.h>
#endif

#ifdef HAVE_MATH_H
 #include<math.h>
#endif

#ifdef HAVE_STDIO_H
 #include<stdio.h>
#endif

#ifdef HAVE_TIME_H
 #include<time.h>
#endif

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

#ifdef WITH_DMALLOC
  #include <dmalloc.h>
#endif

#include"tsnnls.h"

double *random_matrix(int m, int n)

/* Creates a random sparse-ish matrix for test purposes. In each column, we fill in 
   about 1/10 of the entries with nonzero elements. */

{
  double *A = calloc(m*n,sizeof(double));
  int i,j;

  for(j=0;j<n;j++) {

    for(i=0;i<ceil((double)(m)/10.0);i++) {

      A[(rand() % Msize) + Msize*j] = 2.0*((double)(rand())/(double)(RAND_MAX)) - 1.0;

    }

  }

  return A;

}


int main() 
{

  int test,i;
  double *A,*ATT;
  taucs_ccs_matrix *Accs,*ATccs,*ATTccs;
  FILE *outfile;
  int seed;

  seed = 1199383616;
  srand(seed);  

  printf("transpose_test\n"
	 "Random seed is %d.\n"
	 "Running %d tests with %d x %d random matrices...\n\n",
	 seed,NUM_TESTS,Msize,Nsize);

  for(test=0;test<NUM_TESTS;test++) {

    A = random_matrix(Msize,Nsize);
    Accs = taucs_construct_sorted_ccs_matrix(A,Nsize,Msize);

    if (Msize < 10 && Nsize < 10) {

      printf("Test %d.\n"
	     "-----------------------------------------\n"
	     "A = \n",test);
      
      taucs_print_ccs_matrix(Accs);

    }

    ATccs = taucs_ccs_transpose(Accs);

    if (Msize < 10 && Nsize < 10) {

      printf("-----------------------------------------\n"
	     "A^T = \n");
      taucs_print_ccs_matrix(ATccs);

    }
      
    ATTccs = taucs_ccs_transpose(ATccs);

    if (Msize < 10 && Nsize < 10) {

      printf("-----------------------------------------\n"
	     "Comparing A and A^TT....");

    } else {

      printf("%3d. Comparing A and A^TT...",test);

    }

    ATT = taucs_convert_ccs_to_doubles( ATTccs );

    for(i=0;i<Accs->n*Accs->m;i++) {

      if (fabs(ATT[i] - A[i]) > 1e-12) {
	
	printf("FAIL.\n");

	printf("Storing bad matrix A in tprob_A.sparse.\n");
	outfile = fopen("tprob_A.sparse","w");
	taucs_ccs_write_sparse(outfile,Accs);
	fclose(outfile);

	printf("Storing bad matrix AT in tprob_AT.sparse.\n");
	outfile = fopen("tprob_AT.sparse","w");
	taucs_ccs_write_sparse(outfile,ATccs);
	fclose(outfile);

	printf("Storing bad matrix ATT in tprob_ATT.sparse.\n");
	outfile = fopen("tprob_ATT.sparse","w");
	taucs_ccs_write_sparse(outfile,ATTccs);
	fclose(outfile);

        exit(1);

      }

    }

    printf("pass.\n");

    free(A); free(ATT); 
    taucs_ccs_free(Accs);
    taucs_ccs_free(ATccs);
    taucs_ccs_free(ATTccs);
    
  }

  exit(0);
}
  

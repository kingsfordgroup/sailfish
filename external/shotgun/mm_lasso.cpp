/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Alogregprob->ny KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.

   Wrapper code for running shotgun lasso using matrix market format.
   Written by Danny Bickson, CMU
*/
#include <stdio.h>

#include <unistd.h>
#include "common.h"
void usage(char * argv[]){
       fprintf(stderr, "Usage: %s\n\t-m matrix A in sparse matrix market format\n"
                       "\t -v vector y in sparse matrix market format\n"
		       "\t -o output file name (will contain solution vector x, default is x.mtx)\n"
	               "\t -a algorithm (1=lasso, 2=logitic regresion, 3 = find min lambda for all zero solution)\n"
		       "\t -t convergence threshold (default 1e-5)\n"
	               "\t -k solution path length (for lasso)\n"
	               "\t -i  max_iter (default 100)\n"
	               "\t -n num_threads (default 2)\n"
		       "\t -l lammbda - positive weight constant (default 1)\n"
		       "\t -V verbose: 1=verbose, 0=quiet (default 0) \n", argv[0]);  
       exit(1);

}


int main(int argc, char * argv[]){

   bool inputmat = false; bool inputvec = false;  
    double threshold = 1e-5;
    int K = 0;
    int verbose = 0;
    int maxiter = 100;
    int numthreads = 2;
    double lambda = 1;
    int algorithm = 1;
     char c;
     char matfile[256]={0};
     char vecfile[256]={0};
     char outfile[256]={0};
     strcpy(outfile, "x.mtx");

     while ((c = getopt (argc, argv, "m:v:o:a:t:k:i:n:l:V:")) != -1)
         switch (c)
           {
           case 'm':
             strncpy(matfile, optarg, 256);
	     inputmat = true;
             break;
           case 'v':
             strncpy(vecfile, optarg, 256);
	     inputvec = true;
             break;
           case 'o':
             strncpy(outfile, optarg, 256);
	     break;
	   case 'a':	
	     algorithm = atoi(optarg);
	     if (algorithm < 1 || algorithm > 3){
	       fprintf(stderr, "Algorithm should be 1-3\n");
	       exit(1);
	     }
	     break;
           case 't':	
	      threshold = atof(optarg);	
              break;
	   case 'k':	 
	      K = atoi(optarg);
	      break;
	   case 'i':
	      maxiter=atoi(optarg);
	      break;
	   case 'n':
	     numthreads = atoi(optarg);
	     if (numthreads < 1){
		fprintf(stderr, "Wrong number of threads %d\n", numthreads);
	        exit(1);
	     }
	     break;
           case 'l':
	     lambda=atof(optarg);
	     break;
	   case 'V':
	     verbose = atoi(optarg);
	     break;
           default:
             usage(argv);
	   }
     
    if (!inputvec || !inputmat){
	fprintf(stderr, "Matrix and vector files are mandaroty inputs\n");
	usage(argv);
    }

    shotgun_data prob;
    convert_2_mat(matfile, &prob);
    
    
    convert_2_vec(vecfile, &prob);

    #ifndef DISABLE_OMP
    omp_set_num_threads(numthreads);
    mexPrintf("OMP threads = %d\n", numthreads);
    #endif
    mexPrintf("d = %d, n=%d, y-dim=%u\n", prob.nx,prob.ny, (unsigned int)prob.y.size());
    if ((int)prob.y.size() != prob.ny){
       fprintf(stderr, "Error: wrong size of observation vector\n");
       exit(1);
    }
   
    srand(time(NULL));


    bool all_zero = false;
    if (algorithm == 1) {
         solveLasso(&prob, lambda, K, threshold, maxiter, verbose);
    } else if (algorithm == 2) {
         compute_logreg(&prob, lambda, threshold, maxiter, verbose, all_zero);
    } else if (algorithm == 3){
	double max_lambda = 1e10;
	double min_lambda = 1e-10;
	int tries = 0;
	lambda = max_lambda;
        while (tries < 50){
	   printf("Trying out lambda %g", lambda);
           compute_logreg(&prob, lambda, threshold, 3, verbose, all_zero);
           if (all_zero){
                max_lambda = lambda;
		lambda = (max_lambda + min_lambda) / 2;
		printf("OK\n");
                all_zero = false;
           }
	   else {
		min_lambda = lambda;
		lambda = (max_lambda + min_lambda) / 2;
		printf("Too small\n");
           }
	   tries++;
        }
       printf("Found minimal lambda with all zero solution %g\n", lambda);
    }
    else {
     usage(argv);
    }
      
    int * I = new int[prob.nx];
    int * J = new int[prob.nx];
    double * val = new double[prob.nx];
    for (int i=0; i< prob.nx; i++){
      I[i] = i;
      J[i] = 0;
      val[i] = prob.x[i];
    } 
   write_to_file(outfile, I,J,val,prob.nx,1,prob.nx);
}

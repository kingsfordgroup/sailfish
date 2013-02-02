/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/


//
// Common datastructures for lasso and logistic regression
// We call values of vector x "features" and "weights"
//

#ifndef __SHOTGUN_COMMON_
#define __SHOTGUN_COMMON_


#include <cmath>
#include <iostream>
#include <memory.h>
#include <cstdio>
#include <string>
#include <vector>
#include <queue>
#include <list>
#include <map>

#ifndef MEX
#define mexPrintf printf
#else
#include "mex.h"
#endif


#ifndef DISABLE_OMP
#include <omp.h>
#endif 

#include "cas_array.h" // Concurrent array

typedef double valuetype_t;

inline bool isnan_(double x) {
    return !(x<=0.0 || x>0);
}

#define sqr(x) ((x)*(x))


// Simple sparse array
struct sparse_array {
    std::vector<unsigned int> idxs;
    std::vector<valuetype_t> values;   
    
    int length() {
        return idxs.size();
    }
    void add(unsigned int _idx, valuetype_t value) {
        idxs.push_back(_idx);
        values.push_back(value);
    }
};



struct feature {
    valuetype_t covar;  // Covariane of a a column
    valuetype_t Ay_i;   // (Ay)_i
};





struct shotgun_data {
    // Column-wise sparse matrix representation. 
    std::vector<sparse_array> A_cols, A_rows;
    
    // Vector of observations
    std::vector<valuetype_t> y;
    std::vector<valuetype_t> x; 
    std::vector<feature> feature_consts;
    
    // Problem size and configuration
    int nx;  // Number of features/weights
    int ny;  // Number of datapoints
    
    // Lookup tables for optimization
    cas_array<valuetype_t> Ax;
    cas_array<valuetype_t> expAx;
    cas_array<double> Gmax;
};

void convert_2_mat(const char * filename, shotgun_data * prob);
void convert_2_vec(const char * filename, shotgun_data * prob);

double solveLasso(shotgun_data * prob, double lambda, int regpathlength, double threshold, int maxiter, int verbose);
void compute_logreg(shotgun_data * prob, double lambda, double term_threshold, int max_iter, int verbose, bool &all_zero);
void write_to_file(const char * filename, int * I, int * J, double * val, int M, int N, int nnz);

#endif


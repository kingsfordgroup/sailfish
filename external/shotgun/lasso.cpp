/*
   Copyright [2011] [Aapo Kyrola, Joseph Bradley, Danny Bickson, Carlos Guestrin / Carnegie Mellon University]

   Licensed under the Apache License, Version 2.0 (the "License");
   you may not use this file except in compliance with the License.
   You may obtain a copy of the License at

       http://www.apache.org/licenses/LICENSE-2.0

   Unless required by applicable law or agreed to in writing, software
   distributed under the License is distributed on an "AS IS" BASIS,
   WITHOUT WARRANTIES OR CONDITIONS OF Alassoprob->ny KIND, either express or implied.
   See the License for the specific language governing permissions and
   limitations under the License.
*/
// Optimization problem
//      \arg \min_x ||Ax-y|| + \lambda |x|_1
//

#include "common.h"

// Problem definition
shotgun_data * lassoprob;


// Major optimization: always keep updated vector 'Ax'
 
void initialize_feature(int feat_idx) {
    sparse_array& col = lassoprob->A_cols[feat_idx];
    feature& feat = lassoprob->feature_consts[feat_idx];
    lassoprob->x[feat_idx] = 0.0;

    // Precompute covariance of a feature
    feat.covar = 0;
    for(int i=0; i<col.length(); i++) {
        feat.covar += sqr(col.values[i]);
    }
    feat.covar *= 2;
    
    // Precompute (Ay)_i
    feat.Ay_i = 0.0;
    for(int i=0; i<col.length(); i++) {
        feat.Ay_i += col.values[i]*lassoprob->y[col.idxs[i]];
    }
    feat.Ay_i *= 2;

}

void initialize() {
    lassoprob->feature_consts.reserve(lassoprob->nx);
    lassoprob->x.resize(lassoprob->nx);
    lassoprob->Ax.resize(lassoprob->ny);
    #pragma omp for
    for(int i=0; i<lassoprob->nx; i++) {
        initialize_feature(i);
    }
}

valuetype_t soft_thresholdO(valuetype_t _lambda, valuetype_t shootDiff) {
    return (shootDiff > _lambda)* (_lambda - shootDiff) + 
	               (shootDiff < -_lambda) * (-_lambda- shootDiff) ;
}

valuetype_t soft_threshold(valuetype_t _lambda, valuetype_t shootDiff) {
  if (shootDiff > _lambda) return _lambda - shootDiff;
  if (shootDiff < -_lambda) return -_lambda - shootDiff ;
  else return 0;
}

double shoot(int x_i, valuetype_t lambda) {
    feature& feat = lassoprob->feature_consts[x_i];
    valuetype_t oldvalue = lassoprob->x[x_i];
    
    if (feat.covar == 0.0) return 0.0; // Zero-column
    
    // Compute dotproduct A'_i*(Ax)
    valuetype_t AtAxj = 0.0;
    sparse_array& col = lassoprob->A_cols[x_i];
    int len=col.length();
    for(int i=0; i<len; i++) {
        AtAxj += col.values[i] * lassoprob->Ax[col.idxs[i]];
    }
    
    valuetype_t S_j = 2 * AtAxj - feat.covar * oldvalue - feat.Ay_i;
    valuetype_t newvalue = soft_threshold(lambda,S_j)/feat.covar;
    valuetype_t delta = (newvalue - oldvalue);
    
    // Update ax
    if (delta != 0.0) {
        for(int i=0; i<len; i++) {
	       lassoprob->Ax.add(col.idxs[i], col.values[i] * delta);
        }
        
        lassoprob->x[x_i] = newvalue;
    }
    return std::abs(delta);
}

// Find such lambda that if used for regularization,
// optimum would have all weights zero.
valuetype_t compute_max_lambda() {
    valuetype_t maxlambda = 0.0;
    for(int i=0; i<lassoprob->nx; i++) {
        maxlambda = std::max(maxlambda, std::abs(lassoprob->feature_consts[i].Ay_i));
    }
    return maxlambda;
}

valuetype_t get_term_threshold(int k, int K, double delta_threshold) {
  // Stricter termination threshold for last step in the optimization.
  return (k == 0 ? delta_threshold  : (delta_threshold + k*(delta_threshold*50)/K));
}

 valuetype_t compute_objective(valuetype_t _lambda, std::vector<valuetype_t>& x, double & l0x, valuetype_t * l1x = NULL, valuetype_t * l2err = NULL) {
    double least_sqr = 0;
     
    for(int i=0; i<lassoprob->ny; i++) {
        least_sqr += (lassoprob->Ax[i]-lassoprob->y[i])*(lassoprob->Ax[i]-lassoprob->y[i]);					
    }
    
     // Penalty check for feature 0
    double penalty = 0.0;
    for(int i=0; i<lassoprob->nx; i++) {
        penalty += std::abs(lassoprob->x[i]);
        l0x += (lassoprob->x[i] == 0);
    }
    if (l1x != NULL) *l1x = penalty;
    if (l2err != NULL) *l2err = least_sqr;
    return penalty * _lambda + least_sqr;
}



void main_optimization_loop(double lambda, int regpathlength, double threshold, int maxiter, int verbose) {
    // Empirically found heuristic for deciding how malassoprob->ny steps
    // to take on the regularization path
    int regularization_path_length = (regpathlength <= 0 ? 1+(int)(lassoprob->nx/2000) : regpathlength);
    mexPrintf("regularization_path_length = %d\n",regularization_path_length);

    valuetype_t lambda_max = compute_max_lambda();
    valuetype_t lambda_min = lambda;
    valuetype_t alpha = pow(lambda_max/lambda_min, 1.0/(1.0*regularization_path_length));
    int regularization_path_step = regularization_path_length;

    double delta_threshold = threshold;
    long long int num_of_shoots = 0;
    int counter = 0;
    int iterations = 0;
    double *delta = new double[lassoprob->nx];
    do {
        ++iterations;
        if (iterations >= maxiter && maxiter > 0) {
            mexPrintf("Exceeded max iterations: %d", maxiter);
            return;
        }
        double maxChange=0;
        lambda = lambda_min * pow(alpha, regularization_path_step);
        
        // Parallel loop
        #pragma omp parallel for  
        for(int i=0; i<lassoprob->nx; i++) 
            delta[i] = shoot(i, lambda);
        maxChange = 0;
        // TODO: use OpenMP reductions (although - not much to gain)
        for(int i=0; i<lassoprob->nx; i++)
            maxChange = (maxChange < delta[i] ? delta[i] : maxChange);
        num_of_shoots += lassoprob->nx;
        counter++;

        // Convergence check.
        // We use a simple trick to converge faster for the intermediate sub-optimization problems
        // on the regularization path. This is important because we do not care about accuracy for
        // the intermediate problems, just want a good warm start for next round.
        bool converged = (maxChange <= get_term_threshold(regularization_path_step,regularization_path_length,delta_threshold));
        if (converged || counter>std::min(100, (100-regularization_path_step)*2)) {
            counter = 0;
            regularization_path_step--; 
        }
        if (verbose){
          double l1x = 0, l2err = 0, l0x = 0;
          valuetype_t obj = compute_objective(lambda, lassoprob->x, l0x, &l1x, &l2err);
          mexPrintf("Objective: %g L1 : %g L2err: %g L0 %d\n", obj, l1x, l2err, (int)l0x);
        }
    } while (regularization_path_step >= 0);
    delete[] delta;
    mexPrintf("Num of shoots = %lld\n", num_of_shoots);
}


double solveLasso(shotgun_data  * probdef, double lambda, int regpathlength, double threshold, int maxiter, int verbose) {
    lassoprob = probdef;
    initialize();
    main_optimization_loop(lambda, regpathlength, threshold, maxiter, verbose); 
    return 0;
}



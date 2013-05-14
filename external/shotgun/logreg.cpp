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
*/
//
// Implementation of Shotgun Logreg - a parallel Logistic Lasso solver - for OpenMP.
// Optimization problem
//      \arg \max_x \sum log(1+exp(-yAx) - \lambda |x|_1
//
// Based on Coordinate Descent Newton algorithm described in
// Yuan, Chang et al. : 
//      A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification
//  
// \author Aapo Kyrola akyrola@cs.cmu.edu
// modified by Danny Bickson, CMU 
 
#include "common.h"


shotgun_data * logregprob;

// Parameters
double cdn_beta = 0.5;
double cdn_sigma = 0.01;

double Gmax;
double Gmax_old;
double Gmax_init;
double *  xjneg;
int pos_y, neg_y;
bool * active;
bool shuffle = true;

int test_ratio = 10;
double compute_objective_logreg(double lambda, double * l1x = NULL, double * loglikelihood = NULL, int * _l0 = NULL,
                double * testobj = NULL);


inline double sign(double a) {
    return (a<0 ? -1 : (a == 0 ? 0 : 1));
}
inline void swap(int &t1, int &t2) { int tmp=t2; t2=t1; t1=tmp; }
  
   
  

// Computes L_j'(0) and L_j''
// See: Yuan, Chang et al. : 
//  A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification
//  equation (25)
inline void logreg_cdn_derivandH(double lambda, int x_i, double &G, double& H) {
    G = H = 0;
    sparse_array& col = logregprob->A_cols[x_i];
    assert(col.length()>0);
    for(int i=0; i<col.length(); i++) {
       int rowi = col.idxs[i];
       double val = col.values[i];
 
       double exp_wTxind = logregprob->expAx[rowi];
	double tmp1 = val/(1+exp_wTxind);
	double tmp2 =  tmp1;
	double tmp3 = tmp2*exp_wTxind;
	G += tmp2;
	H += tmp1*tmp3;
    }
    G = -G +  xjneg[x_i];
    G /= lambda;
    H /= lambda;
    if (H<1e-5) {
        H = 1e-5;
    }
 }
 


 // L_j(x + diff*e_j)-L_j(x) (see function g_xi())
 // (eq. 18)
inline double logreg_cdn_Ldiff(double lambda, int x_i, double diff) {
    double sum = 0.0;
    sparse_array& col = logregprob->A_cols[x_i];
    for(int i=0; i<col.length(); i++) {
       int rowi = col.idxs[i];
       double dc = diff * col.values[i];
       double expdiff = exp(dc);
       double expAXdiff = logregprob->expAx[rowi] * expdiff;
       //assert(!isnan(logregprob->expAx[rowi]));
       //assert(expAXdiff + expdiff != 0);
       double ds = log((expAXdiff+1)/(expAXdiff+expdiff));
       //assert(!isnan(ds));
       sum +=  ds;
    } 
    if (isnan(sum)){
        fprintf(stderr, "Got numerical error: please verify that in your dataset there are no columns of matrix A with all zeros. Encountered error in column %d\n", x_i);
        exit(1);
    }
    return 1.0/lambda * (diff*xjneg[x_i] + sum);
}
 
 // Equation 17. One-variable function that equals the change in loss function
 // when x_i is change by z
double g_xi(double z, int x_i, double lambda) {
    double xv = logregprob->x[x_i];
    return std::abs(xv+z) - std::abs(xv) + logreg_cdn_Ldiff(lambda, x_i, z);
}
 

// Compute A'A in parallel. Note: could also compute in demand?
void initialize_all() {
    logregprob->x.resize(logregprob->nx);
    logregprob->Ax.resize(logregprob->ny);
    logregprob->expAx.resize(logregprob->ny);
    logregprob->Gmax.resize(1);

    active = (bool *) calloc(logregprob->nx,sizeof(bool));
    xjneg = (double *) calloc(logregprob->nx,sizeof(double));
    pos_y = 0;
    for(int i=0; i<logregprob->ny; i++)  {
        logregprob->expAx[i] = 1.0; // since(exp(0) = 1)
        pos_y += (logregprob->y[i] == 1);
    }
    neg_y = logregprob->ny-pos_y;
    
    assert(pos_y > 0 && neg_y > 0);
    for(int i=0; i<logregprob->nx; i++) active[i] = true;
    
    // use the trick that log(1+exp(x))=log(1+exp(-x))+x
    #pragma omp parallel for
    for(int i=0; i<logregprob->nx; i++) {
        sparse_array& col = logregprob->A_cols[i];
        for(int j=0; j<col.length(); j++) {
            if (logregprob->y[col.idxs[j]] == -1)
                xjneg[i] += col.values[j];
        }
    }
}

void recompute_expAx() {
    #pragma omp parallel for
    for(int i=0; i<logregprob->ny; i++) {
        double Ax=0;
        sparse_array &row = logregprob->A_rows[i];
        for(int j=0; j<row.length(); j++) {
            Ax += logregprob->x[row.idxs[j]]*row.values[j];
        }
        logregprob->expAx[i] = exp(Ax);
    }
}
 
 
// Yua, Chang, Hsieh and Lin: A Comparison of Optimization Methods and Software for Large-scale L1-regularized Linear Classification; p. 14
double shoot_cdn(int x_i, double lambda) {
    // Compute d: (equation 29), i.e the solution to the quadratic approximation of the function 
    // for weight x_i
    if (!active[x_i]) return 0.0;
    
    double violation = 0.0;
    double xv = logregprob->x[x_i];
    
    double Ld0, Ldd0;
    logreg_cdn_derivandH(lambda, x_i, Ld0, Ldd0);    
    
    double Gp = (Ld0+1);
    double Gn = (Ld0-1);

    // Shrinking (practically copied from LibLinear)
    if (xv == 0) {
        if (Gp<0) {
            violation = -Gp;   
        } else if (Gn > 0) {
            violation = Gn;
       	} else if(Gp>Gmax_old/logregprob->ny && Gn<-Gmax_old/logregprob->ny) {
            // Remove
            active[x_i] = false;
            return 0.0;
        }   
    } else if(xv > 0)
      violation = fabs(Gp);
    else
      violation = fabs(Gn);
    
    // TODO: should use a reduction here! Or lock.
    //if (Gmax < violation)
    //    Gmax = violation;
    logregprob->Gmax.max(0, violation);
    //printf("node %d violation %g Ld0 %g Ldd0 %g Gp %g Gn %g xv %g \n", x_i, violation, Ld0, Ldd0, Gp, Gn, xv);
    // Newton direction d
    double rhs = Ldd0*xv;
    double d;
    if (Gp<= rhs) {
        d = -(Gp/Ldd0);
    } else if (Gn >= rhs) {
        d = -(Gn/Ldd0); 
    } else {
        d = -xv;
    }
    
     if (std::abs(d) < 1e-12) {
        return 0.0;
    }
    // Small optimization
    d = std::min(std::max(d,-10.0),10.0);
   
    // Backtracking line search (with Armijo-type of condition) (eq. 30)
    int iter=0;
    int max_num_linesearch=20;
    double gamma = 1.0; // Note: Yuan, Chang et al. use lambda instead of gamma
    double delta = (Ld0 * d + std::abs(xv+d) - std::abs(xv));
    double rhs_c = cdn_sigma * delta;    
    do {
        double change_in_obj = g_xi(d, x_i, lambda);
        if (change_in_obj <= gamma * rhs_c) {
            // Found ok.
            logregprob->x[x_i] += d;
            // Update dot products (Ax)
            sparse_array &col = logregprob->A_cols[x_i];
            #pragma omp parallel for
            for(int i=0; i<col.length(); i++) {
                logregprob->expAx.mul(col.idxs[i], exp(d * col.values[i]));
            }
            return std::abs(d);
        }
        gamma *= 0.5;
        d *= 0.5;
    } while(++iter < max_num_linesearch); 
    recompute_expAx();
    return 0.0;
}

 
  
 

 
/** 
  * Execute "shoot" update of a feature.
  */
void shoot_logreg(int x_i, double lambda) {
    // Some columns may be empty:
    if (logregprob->A_cols[x_i].length() == 0) return;
    shoot_cdn(x_i, lambda);
}


/**
  * Main optimization loop.
  * Note: this version of code does NOT have special version for sequential
  * version. Sequential version is slightly lighter than parallel version because
  * it can use a more efficient operation to maintain the active set. In practice, this
  * has little effect. Also, it does not need to have a atomic array for maintaining Ax.
  * For the experiments in the paper, a special sequential code was used for fairness.
  */
void compute_logreg(shotgun_data * prob, double lambda, double term_threshold, int max_iter, int verbose, bool & all_zero) {
    all_zero = false;
    logregprob = prob;
    //double l1x, loglikelihood;
    int iterations = 0;//, t=0;
    long long int num_of_shoots = 0;

    std::vector<int> shuffled_indices;
    for(int j=0; j<logregprob->nx; j++) shuffled_indices.push_back(j);

    Gmax_old = 1e30;
    // Adjust threshold similarly as liblinear
    initialize_all();
    term_threshold =  term_threshold*std::min(pos_y,neg_y)/double(logregprob->ny);
    
    while(true) {
        int active_size = logregprob->nx;
        num_of_shoots += active_size;   
        
         // Randomization
        if (shuffle)
        for(int j=0; j<active_size; j++) {
                int i = j+rand()%(active_size-j);
                swap(shuffled_indices[i], shuffled_indices[j]);
        } 
            
            /* Main parallel loop */
        #pragma omp parallel for
        for(int s=0; s<active_size; s++) {
            int x_i = shuffled_indices[s];
            shoot_logreg(x_i, lambda);
        }
            
        /* Gmax handling */
        Gmax_old = logregprob->Gmax[0];
        if (iterations == 0) {
            Gmax_init = Gmax_old;
        }
        
        iterations++;
        
        //std::cout << Gmax.get_value() <<  " " << Gmax_init << " " <<  term_threshold*Gmax_init << std::endl;
        if (iterations > max_iter && max_iter>0) {
            mexPrintf("Exceeded max iterations: %d\n", max_iter);
            break;
        }
        
        for(int i=0; i<logregprob->nx; i++) shuffled_indices[i] = i;
        active_size = logregprob->nx;
        for(int s=0; s<active_size; s++) {
            int j = shuffled_indices[s];
            if (!active[j]) {
                active_size--;
                swap(shuffled_indices[s], shuffled_indices[active_size]);
                s--;
            }
        }
        
        if (logregprob->Gmax[0] <= term_threshold*Gmax_init) {
           // std::cout << active_size << std::endl;
            if (active_size == logregprob->nx) {
                printf("Encountered all zero solution! try to decrease lambda\n");
 		all_zero = true;
   		break;              
            } else {
                 Gmax_old = 1e30;
                 for(int i=0; i<logregprob->nx; i++) active[i] = true;
                 active_size=logregprob->nx;
                 recompute_expAx();
                 //continue;
            }
        }  

     if (verbose){
         double l1x=0, loglikelihood=0;
         int l0=0;
         double obj = compute_objective_logreg(lambda, &l1x, &loglikelihood, &l0, NULL);
        printf("objective is: %g l1: %g loglikelihood %g l0: %d\n", obj, l1x, loglikelihood, l0); 
	if (l1x == 0)	
	    all_zero = true;
     }
   }// end iterations

   if (!verbose){
      double l1x=0, loglikelihood=0;
      int l0=0;
      double obj = compute_objective_logreg(lambda, &l1x, &loglikelihood, &l0, NULL);
      if (l1x == 0)	
	all_zero = true;
      printf("objective is: %g l1: %g loglikelihood %g l0: %d\n", obj, l1x, loglikelihood, l0); 
   }

   delete[] active;
   delete[] xjneg;
  mexPrintf("Finished Shotgun CDN in %d iterations\n", iterations);
}




#include "common.h"
extern shotgun_data * logregprob;
double logloss(double x) {
     if (x > (-10) && x < 10) {
        return log(1 + exp(x));
     } else if (x<= (-10)) {
        return 0.0;
     } else {
        return x;
     }
 }

double compute_llhood(std::vector<double>& xx, std::vector<double> & yy, std::vector<sparse_array> & rows) {
     double llhood = 0;
     int rs = logregprob->ny;
     for(int i=0; i<rs; i++) {
        sparse_array& row = rows[i];
        double Ax=0;
        for(int j=0; j<row.length(); j++) {
            Ax += xx[row.idxs[j]]*row.values[j];
        }
        llhood -= logloss(-yy[i]*Ax);
     }
     return llhood;
 }

 double compute_objective_logreg(double lambda, double * l1x = NULL, double * loglikelihood = NULL, int * _l0 = NULL,
                double * testobj = NULL) {
    std::vector<double> & x = logregprob->x;
    double penalty = 0.0;
    int l0 = 0;
    for(int j=0; j<logregprob->nx; j++) {
        penalty += std::abs(x[j]);
        l0 += (x[j] != 0);
    }
    double llhood = compute_llhood(logregprob->x, logregprob->y, logregprob->A_rows);
    double llhood_test = 0;// //TODO compute_llhood(logregprob->x, logregprob->ytest, logregprob->Atest_rows);
    if (l1x != NULL) *l1x = penalty;
    if (loglikelihood != NULL) *loglikelihood = llhood;
    if (_l0 != NULL) *_l0 = l0;
    if (testobj != NULL) *testobj = -penalty*lambda + llhood_test; 
    return -penalty*lambda + llhood;
}


